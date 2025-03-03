using LinearAlgebra
using SparseArrays
using GLMakie

struct Material
    E::Float64
    ν::Float64
    thickness::Float64
end

struct CohesiveProperties
    K0::Float64      # Cohesive stiffness for displacement (N/m³)
    K_rot::Float64   # Cohesive stiffness for rotation (N·m/m)
    σ_max::Float64   # Maximum traction stress (Pa)
    δ_c::Float64     # Critical separation (m)
    m::Float64       # Fatigue degradation rate
end

const Si = Material(170e9, 0.28, 675e-6)
const Parylene = Material(2.8e9, 0.4, 10e-6)
const Steel = Material(200e9, 0.3, 200e-6)

const CZM = CohesiveProperties(1e8, 1e6, 5e6, 5e-6, 1e-6)

const L_steel = 60e-3
const L_layered = 2e-3
const n_elem = 100
const dx = L_steel / n_elem
const nodes = collect(0:dx:L_steel)

const left_start_pos = 20e-3
const left_end_pos = left_start_pos + L_layered
const left_start_elem = floor(Int, left_start_pos / dx) + 1
const left_end_elem = floor(Int, left_end_pos / dx)
const n_elem_left = left_end_elem - left_start_elem + 1

function beam_element(E, I, L)
    k = E * I / L^3 * [
        12   6*L   -12   6*L
        6*L  4*L^2 -6*L  2*L^2
        -12  -6*L   12  -6*L
        6*L  2*L^2 -6*L  4*L^2
    ]
    return k
end

function assemble_steel_system(force::Float64)
    n_dof = length(nodes) * 2
    K = spzeros(n_dof, n_dof)
    F = zeros(n_dof)
    
    I_steel = Steel.thickness^3 / 12
    
    for i in 1:n_elem
        idx = [(i-1)*2+1:(i-1)*2+2; (i)*2+1:(i)*2+2]
        k_steel = beam_element(Steel.E, I_steel, dx)
        K[idx, idx] += k_steel
    end
    
    mid = div(n_elem, 2) + 1
    F[mid*2-1] = -force
    
    K[1,1] += 1e8
    K[n_dof-1,n_dof-1] += 1e8
    
    return K, F
end

function assemble_si_system(n_elem_layered::Int)
    n_dof_si = (n_elem_layered + 1) * 2
    K_si = spzeros(n_dof_si, n_dof_si)
    F_si = zeros(n_dof_si)
    
    I_si = Si.thickness^3 / 12
    
    for i in 1:n_elem_layered
        idx_si = [(i-1)*2+1:(i-1)*2+2; i*2+1:i*2+2]
        k_si = beam_element(Si.E, I_si, dx)
        K_si[idx_si, idx_si] += k_si
    end
    
    # Pin both ends for stability
    K_si[1,1] += 1e10
    K_si[end-1,end-1] += 1e10
    
    return K_si, F_si
end

function simulate_fatigue(max_force::Float64, n_cycles::Int)
    damage_left = zeros(n_elem_left + 1)
    separated_left = falses(n_elem_left + 1)
    u_history = Vector{Vector{Float64}}()
    n_dof_steel = length(nodes) * 2
    n_dof_si_left = (n_elem_left + 1) * 2
    n_dof = n_dof_steel + n_dof_si_left
    
    K_steel, F_steel = assemble_steel_system(max_force)
    u_steel_initial = K_steel \ F_steel
    println("Initial Max Steel Deflection (μm): ", maximum(abs.(u_steel_initial[1:2:n_dof_steel])) * 1e6)
    
    K_si_left, F_si_left_base = assemble_si_system(n_elem_left)
    u_initial = zeros(n_dof)
    u_initial[1:n_dof_steel] = u_steel_initial
    
    # Set Si initial deflection and rotation to match steel/Parylene
    offset_left = n_dof_steel
    for i in 1:n_elem_left + 1
        idx_steel_v = (left_start_elem + i - 2) * 2 - 1
        idx_steel_θ = (left_start_elem + i - 2) * 2
        idx_si_v = offset_left + (i-1) * 2 + 1
        idx_si_θ = offset_left + (i-1) * 2 + 2
        u_initial[idx_si_v] = u_steel_initial[idx_steel_v]
        u_initial[idx_si_θ] = u_steel_initial[idx_steel_θ]
    end
    println("Initial Max Si Left Deflection (μm): ", maximum(abs.(u_initial[n_dof_steel+1:2:end])) * 1e6)
    push!(u_history, u_initial)
    
    for cycle in 1:n_cycles
        u_prev = u_history[end]
        u_new = copy(u_prev)
        
        # Steel deflection
        u_new[1:n_dof_steel] = K_steel \ F_steel
        
        # Si left deflection with cohesive forces and moments
        F_si_left = copy(F_si_left_base)
        for i in 1:n_elem_left + 1
            idx_steel_v = (left_start_elem + i - 2) * 2 - 1
            idx_steel_θ = (left_start_elem + i - 2) * 2
            idx_si_v = offset_left + (i-1) * 2 + 1
            idx_si_θ = offset_left + (i-1) * 2 + 2
            
            v_parylene = u_new[idx_steel_v]
            θ_parylene = u_new[idx_steel_θ]
            v_si_prev = u_prev[idx_si_v]
            θ_si_prev = u_prev[idx_si_θ]
            
            δ_v = abs(v_si_prev - v_parylene)
            T = CZM.K0 * δ_v
            damage_left[i] += CZM.m * (δ_v / CZM.δ_c)
            damage_left[i] = min(damage_left[i], 1.0)
            
            if !separated_left[i]
                if T > CZM.σ_max
                    separated_left[i] = true
                    println("Cycle $cycle: Left Node $i separated (T = $T Pa, δ_v = $δ_v m)")
                else
                    F_si_left[(i-1)*2 + 1] = -CZM.K0 * (v_si_prev - v_parylene) * dx
                    F_si_left[(i-1)*2 + 2] = -CZM.K_rot * (θ_si_prev - θ_parylene) * dx
                end
            else
                F_si_left[(i-1)*2 + 1] = -CZM.K0 * δ_v * dx * 0.1
                F_si_left[(i-1)*2 + 2] = 0
            end
        end
        
        u_new[offset_left+1:end] = K_si_left \ F_si_left
        
        push!(u_history, u_new)
    end
    
    return u_history, damage_left
end

function plot_results(u_hist, damage_left)
    n_dof_steel = length(nodes) * 2
    u_last = u_hist[end]
    
    x_steel = nodes .* 1e3
    steel_deflection = u_last[1:2:n_dof_steel] .* 1e6
    
    x_left = range(left_start_pos, left_end_pos, length=n_elem_left + 1) .* 1e3
    
    parylene_left = zeros(n_elem_left + 1)
    si_left = zeros(n_elem_left + 1)
    forces_left = zeros(n_elem_left + 1)
    for i in 1:n_elem_left + 1
        elem = left_start_elem + i - 1
        parylene_left[i] = u_last[(elem-1)*2 + 1] * 1e6
        si_left[i] = u_last[n_dof_steel + (i-1)*2 + 1] * 1e6
        δ = abs(u_last[n_dof_steel + (i-1)*2 + 1] - u_last[(elem-1)*2 + 1])
        forces_left[i] = CZM.K0 * δ * dx
    end
    
    x_separation_left = range(left_start_pos, left_end_pos, length=n_elem_left + 1) .* 1e3
    separation_left = abs.(si_left .- parylene_left)
    
    mid_elem_left = div(n_elem_left, 2) + 1
    mid_idx_steel_left = (left_start_elem + mid_elem_left - 2) * 2 + 1
    mid_idx_si_left = n_dof_steel + (mid_elem_left - 1) * 2 + 1
    crack_left = [abs(u[mid_idx_si_left] - u[mid_idx_steel_left]) * 1e6 for u in u_hist]
    
    cycles = 1:length(u_hist)
    
    println("Max separation left (μm): ", maximum(separation_left))
    println("Max damage left: ", maximum(damage_left))
    println("Max steel deflection (μm): ", maximum(abs.(steel_deflection)))
    println("Max internal force left (N): ", maximum(abs.(forces_left)))
    println("\nSteel Deflection Sample (μm):")
    for i in [1, 25, 34, 50, 75, 101]
        println("Position: ", round(x_steel[i], digits=1), " mm, Deflection: ", round(steel_deflection[i], digits=2))
    end
    println("\nParylene Left Sample (μm):")
    for i in 1:2:length(parylene_left)
        println("Position: ", round(x_left[i], digits=2), " mm, Deflection: ", round(parylene_left[i], digits=2))
    end
    println("\nSi Left Sample (μm):")
    for i in 1:2:length(si_left)
        println("Position: ", round(x_left[i], digits=2), " mm, Deflection: ", round(si_left[i], digits=2))
    end
    println("\nCrack Opening Left Sample (μm, every 100 cycles):")
    for i in 1:100:length(crack_left)
        println("Cycle: ", cycles[i], ", Separation: ", round(crack_left[i], digits=2))
    end
    
    fig = Figure(resolution=(1200, 800))
    
    ax1 = Axis(fig[1, 1], title="Steel Beam Deformation (3-Point Bending)", xlabel="Position (mm)", ylabel="Deflection (μm)")
    lines!(ax1, x_steel, steel_deflection, color=:green, linewidth=2, label="Steel")
    vlines!(ax1, [left_start_pos * 1e3, left_end_pos * 1e3], color=:black, linestyle=:dash, label="Layered Region")
    axislegend(ax1)
    
    ax2 = Axis(fig[2, 1], title="Left Beam Deformation (20–22 mm)", xlabel="Position (mm)", ylabel="Deflection (μm)")
    lines!(ax2, x_left, parylene_left, color=:orange, linewidth=2, label="Parylene (Steel)")
    lines!(ax2, x_left, si_left, color=:purple, linewidth=2, label="Silicon")
    axislegend(ax2)
    
    ax3 = Axis(fig[1, 2], title="Damage Distribution (Left)", xlabel="Position (mm)", ylabel="Damage")
    lines!(ax3, x_separation_left, damage_left, color=:red, linewidth=2)
    
    ax4 = Axis(fig[2, 2], title="Crack Opening vs. Cycles (Left, 21 mm)", xlabel="Cycle", ylabel="Separation (μm)")
    lines!(ax4, cycles, crack_left, color=:blue, linewidth=2)
    
    ax5 = Axis(fig[3, 1:2], title="Internal Forces Between Si and Parylene (Left)", xlabel="Position (mm)", ylabel="Force (N)")
    lines!(ax5, x_separation_left, forces_left, color=:cyan, linewidth=2)
    
    fig[1:3, 1:2] = [ax1 ax3; ax2 ax4; ax5 GridLayout()]
    display(fig)
end

function main()
    max_force = 10.0
    cycles = 1000
    
    println("Starting simulation...")
    u_hist, damage_left = simulate_fatigue(max_force, cycles)
    println("Simulation completed. Generating plots with Makie...")
    
    plot_results(u_hist, damage_left)
    
    return u_hist, damage_left
end

u_hist, damage_left = main()