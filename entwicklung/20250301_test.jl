
###############################

using LinearAlgebra
using SparseArrays
using GLMakie

# Material properties
struct Material
    E::Float64      # Young's modulus (Pa)
    ν::Float64      # Poisson's ratio
    thickness::Float64  # Thickness (m)
end

# Cohesive zone properties
struct CohesiveProperties
    K0::Float64     # Initial stiffness (N/m³)
    σ_max::Float64  # Maximum traction (Pa)
    δ_c::Float64    # Critical separation (m)
    m::Float64      # Fatigue degradation parameter
end

# Define materials
const Si = Material(3e9, 0.28, 675e-6)
const Parylene = Material(2.8e9, 0.4, 10e-6)
const Steel = Material(200e9, 0.3, 200e-6)

const CZM = CohesiveProperties(1e8, 50e6, 5e-6, 0.000001)

const L_steel = 50e-3
const L_layered = 5e-3
const n_elem = 100
const dx = L_steel / n_elem
const nodes = collect(0:dx:L_steel)

const start_pos = 15e-3
const end_pos = start_pos + L_layered
const start_elem = floor(Int, start_pos / dx) + 1
const end_elem = floor(Int, end_pos / dx)
const n_elem_layered = end_elem - start_elem + 1

function beam_element(E, I, L)
    k = E * I / L^3 * [
        12   6*L   -12   6*L
        6*L  4*L^2 -6*L  2*L^2
        -12  -6*L   12  -6*L
        6*L  2*L^2 -6*L  4*L^2
    ]
    return k
end

function cohesive_element(K0, L)
    k = K0 * L / 2 * [
        1  0  -1  0
        0  0   0  0
       -1  0   1  0
        0  0   0  0
    ]
    return k
end

function assemble_system(force::Float64)
    n_dof_steel = length(nodes) * 2
    n_dof_si = (n_elem_layered + 1) * 2
    n_dof = n_dof_steel + n_dof_si
    K = spzeros(n_dof, n_dof)
    F = zeros(n_dof)
    
    I_si = Si.thickness^3 / 12
    I_par = Parylene.thickness^3 / 12
    I_steel = Steel.thickness^3 / 12
    
    for i in 1:n_elem
        idx_steel = [(i-1)*2+1:(i-1)*2+2; (i)*2+1:(i)*2+2]
        k_steel = beam_element(Steel.E, I_steel, dx)
        K[idx_steel, idx_steel] += k_steel
        
        if i >= start_elem && i <= end_elem
            k_par = beam_element(Parylene.E, I_par, dx)
            K[idx_steel, idx_steel] += k_par
        end
    end
    
    for i in start_elem:end_elem
        idx_steel = [(i-1)*2+1:(i-1)*2+2; (i)*2+1:(i)*2+2]
        idx_si = [n_dof_steel + (i-start_elem)*2+1:n_dof_steel + (i-start_elem)*2+2;
                  n_dof_steel + (i-start_elem+1)*2+1:n_dof_steel + (i-start_elem+1)*2+2]
        
        k_si = beam_element(Si.E, I_si, dx)
        K[idx_si, idx_si] += k_si
        
        k_coh = cohesive_element(CZM.K0, dx)
        K[idx_si, idx_si] += k_coh
        K[idx_steel, idx_steel] += k_coh
        K[idx_si, idx_steel] -= k_coh
        K[idx_steel, idx_si] -= k_coh
    end
    
    mid = div(n_elem, 2) + 1
    F[mid*2-1] = -force
    
    K[1,1] += 1e8  # Left steel support
    K[n_dof_steel-1,n_dof_steel-1] += 1e8  # Right steel support
    # Fix Si rotational DOFs at boundaries of layered region to eliminate rank deficiency
    K[n_dof_steel+2,n_dof_steel+2] += 1e8  # Left Si rotation (29 mm)
    K[n_dof-1,n_dof-1] += 1e8  # Right Si rotation (31 mm)
    
    return K, F
end


function fatigue_degradation!(K, u, cycles, damage, n_dof_steel)
    for i in 1:n_elem_layered
        elem = start_elem + i - 1
        idx_steel = (elem-1)*2 + 1
        idx_si = n_dof_steel + (i-1)*2 + 1
        δ = abs(u[idx_si] - u[idx_steel])
        
        if δ > 0
            damage[i] += CZM.m * (δ / CZM.δ_c) * cycles
            damage[i] = min(damage[i], 1.0)
            
            K_factor = (1 - damage[i])
            idx_steel_full = [(elem-1)*2+1:(elem-1)*2+2; (elem)*2+1:(elem)*2+2]
            idx_si_full = [n_dof_steel + (i-1)*2+1:n_dof_steel + (i-1)*2+2;
                          n_dof_steel + i*2+1:n_dof_steel + i*2+2]
            
            k_coh_old = cohesive_element(CZM.K0, dx)
            k_coh_new = cohesive_element(CZM.K0 * K_factor, dx)
            K[idx_si_full, idx_si_full] += (k_coh_new - k_coh_old)
            K[idx_steel_full, idx_steel_full] += (k_coh_new - k_coh_old)
            K[idx_si_full, idx_steel_full] -= (k_coh_new - k_coh_old)
            K[idx_steel_full, idx_si_full] -= (k_coh_new - k_coh_old)
        end
    end
    return K
end

function simulate_fatigue(max_force::Float64, n_cycles::Int)
    damage = zeros(n_elem_layered)
    u_history = Vector{Vector{Float64}}()
    n_dof_steel = length(nodes) * 2  # Define locally
    
    K_base, F = assemble_system(max_force)
    println("Rank of K: ", rank(K_base), " Size: ", size(K_base, 1))
    K = copy(K_base)
    u_initial = K \ F
    println("Initial Max Steel Deflection (μm): ", maximum(abs.(u_initial[1:2:n_dof_steel])) * 1e6)
    println("Initial Max Si Deflection (μm): ", maximum(abs.(u_initial[n_dof_steel+1:2:end])) * 1e6)
    
    for cycle in 1:n_cycles
        u = K \ F
        push!(u_history, copy(u))
        
        K = copy(K_base)
        K = fatigue_degradation!(K, u, 1.0, damage, n_dof_steel)
    end
    
    return u_history, damage
end

function plot_results(u_hist, damage)
    forces=[]
    n_dof_steel = length(nodes) * 2
    u_last = u_hist[end]
    
    x_steel = nodes .* 1e3
    steel_deflection = u_last[1:2:n_dof_steel] .* 1e6
    
    x_layered = range(start_pos, end_pos, length=n_elem_layered + 1) .* 1e3
    parylene_deflection = zeros(n_elem_layered + 1)
    si_deflection = zeros(n_elem_layered + 1)
    for i in 1:n_elem_layered + 1
        elem = start_elem + i - 1
        parylene_deflection[i] = u_last[(elem-1)*2 + 1] * 1e6 - steel_deflection[start_elem]
        si_deflection[i] = u_last[n_dof_steel + (i-1)*2 + 1] * 1e6 

        if i <= n_elem_layered
            δ = abs(u_last[n_dof_steel + (i-1)*2 + 1] - u_last[(elem-1)*2 + 1])
            push!(forces,CZM.K0 * (1 - damage[i]) * δ * dx)  # Force in N (unit width)
        end

    end
    
    x_separation = range(start_pos, end_pos, length=n_elem_layered) .* 1e3
    separation = zeros(n_elem_layered)
    for i in 1:n_elem_layered
        elem = start_elem + i - 1
        separation[i] = abs(u_last[n_dof_steel + (i-1)*2 + 1] - u_last[(elem-1)*2 + 1]) * 1e6
    end
    
    mid_elem = div(n_elem_layered, 2) + 1
    mid_idx_steel = (start_elem + mid_elem - 2) * 2 + 1
    mid_idx_si = n_dof_steel + (mid_elem - 1) * 2 + 1
    crack_over_cycles = [abs(u[mid_idx_si] - u[mid_idx_steel]) * 1e6 for u in u_hist]
    cycles = 1:length(u_hist)
    
    println("Max separation (μm): ", maximum(separation))
    println("Max damage: ", maximum(damage))
    println("Max steel deflection (μm): ", maximum(abs.(steel_deflection)))
    println("\nSteel Deflection Sample (μm):")
    for i in [1, 25, 50, 75, 101]
        println("Position: ", round(x_steel[i], digits=1), " mm, Deflection: ", round(steel_deflection[i], digits=2))
    end
    println("\nParylene Deflection Sample (μm):")
    for i in 1:5:length(parylene_deflection)
        println("Position: ", round(x_layered[i], digits=2), " mm, Deflection: ", round(parylene_deflection[i], digits=2))
    end
    println("\nSi Deflection Sample (μm):")
    for i in 1:5:length(si_deflection)
        println("Position: ", round(x_layered[i], digits=2), " mm, Deflection: ", round(si_deflection[i], digits=2))
    end
    println("\nCrack Opening Sample (μm, every 100 cycles):")
    for i in 1:100:length(crack_over_cycles)
        println("Cycle: ", cycles[i], ", Separation: ", round(crack_over_cycles[i], digits=2))
    end
    
    fig = Figure(resolution=(1200, 800))
    
    fig = Figure(resolution=(1200, 1000))  # Increased height for extra plot
    
    ax1 = Axis(fig[1, 1], title="Steel Beam Deformation (3-Point Bending)", xlabel="Position (mm)", ylabel="Deflection (μm)")
    lines!(ax1, x_steel, steel_deflection, color=:green, linewidth=2, label="Steel")
    vlines!(ax1, [start_pos * 1e3, end_pos * 1e3], color=:black, linestyle=:dash, label="Layered Region")
    axislegend(ax1)
    
    ax2 = Axis(fig[2, 1], title="Two-Layered Beam Deformation (29–31 mm)", xlabel="Position (mm)", ylabel="Deflection (μm)")
    lines!(ax2, x_layered, parylene_deflection, color=:orange, linewidth=2, label="Parylene (Steel)")
    lines!(ax2, x_layered, si_deflection, color=:purple, linewidth=2, label="Silicon")
    axislegend(ax2)
    
    ax3 = Axis(fig[1, 2], title="Damage Distribution", xlabel="Position (mm)", ylabel="Damage")
    lines!(ax3, x_separation, damage, color=:red, linewidth=2)
    
    ax4 = Axis(fig[2, 2], title="Crack Opening at Midpoint (30 mm) vs. Cycles", xlabel="Cycle", ylabel="Separation (μm)")
    lines!(ax4, cycles, crack_over_cycles, color=:blue, linewidth=2)
    
    ax5 = Axis(fig[3, 1:2], title="Internal Forces Between Si and Parylene", xlabel="Position (mm)", ylabel="Force (N)")
    lines!(ax5, x_separation, forces, color=:cyan, linewidth=2)
    
    fig[1:3, 1:2] = [ax1 ax3; ax2 ax4; ax5 GridLayout()]
    display(fig)
end

function main()
    max_force = -10.0
    cycles = 10000
    
    println("Starting simulation...")
    u_hist, damage = simulate_fatigue(max_force, cycles)
    println("Simulation completed. Generating plots with Makie...")
    
    plot_results(u_hist, damage)
    
    return u_hist, damage
end

u_hist, damage = main()

