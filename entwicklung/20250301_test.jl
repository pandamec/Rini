
###############################
using Rini
using GLMakie
# Eigenschaften der Materialien
struct Material
    E::Float64          # Young's modulus (Pa)
    ν::Float64          #   Poisson's ratio
    thickness::Float64  # Thickness (m)
end

# CZM Eigenschaften
struct CohesiveProperties
    K0::Float64     # Initial stiffness (N/m³)
    σ_max::Float64  # Maximum traction (Pa)
    δ_c::Float64    # Critical separation (m)
    m::Float64      # Fatigue degradation parameter
end

struct TestSetup
    L_steel::Float64  
    L_layered::Float64  
    n_elem::Int 
    dx::Float64  
    nodes::Vector{Float64}
    start_pos::Float64  
    end_pos::Float64  
    start_elem::Int
    end_elem::Int
    n_elem_layered::Int
end

# Materialien der Layers
const Si        =   Material(3e9, 0.28, 675e-6)
const Parylene  =   Material(2.8e9, 0.4, 10e-6)
const Steel     =   Material(200e9, 0.3, 200e-6)
const CZM       =   CohesiveProperties(1e8, 30e6, 0.0035, 0.000001)

# Geometrie der Probe
const L_steel   = 60e-3
const L_layered = 5e-3
const n_elem    = 200               # Finite Elemente Methode benutzt fuer die Modellierung
const dx        = L_steel / n_elem
const nodes     = collect(0:dx:L_steel)

const start_pos         = 15e-3
const end_pos           = start_pos + L_layered
const start_elem        = floor(Int, start_pos / dx) + 1
const end_elem          = floor(Int, end_pos / dx)
const n_elem_layered    = end_elem - start_elem + 1

const setup             = TestSetup(L_steel,L_layered,n_elem,dx,nodes,start_pos,end_pos,start_elem,end_elem,n_elem_layered)
# Test setup
max_force = -10.0
cycles = 10
    
println("Starting simulation...")
u_hist, damage = simulate_fatigue(setup,max_force, cycles,Si,Parylene,Steel,CZM)
println("Simulation completed. Generating plots with Makie...")

function plot_results(u_hist, damage)
    forces=[]
    n_dof_steel = length(nodes) * 2
    u_last = u_hist[end]
    
    x_steel = nodes .* 1e3
    steel_deflection = u_last[1:2:n_dof_steel] .* 1e6
    
    x_layered = range(start_pos, end_pos, length=n_elem_layered + 1) .* 1e3
    parylene_deflection = zeros(n_elem_layered + 1)
    si_deflection = zeros(n_elem_layered + 1)

    ## Interface Krafte
    for i in 1:n_elem_layered + 1
        elem = start_elem + i - 1
        parylene_deflection[i] = u_last[(elem-1)*2 + 1] * 1e6 #- steel_deflection[start_elem]
        si_deflection[i] = u_last[n_dof_steel + (i-1)*2 + 1] * 1e6 + steel_deflection[elem]

        if i <= n_elem_layered
            δ = abs(u_last[n_dof_steel + (i-1)*2 + 1] - u_last[(elem-1)*2 + 1])
            
            push!(forces,CZM.K0 * (1 - damage[i]) * δ * dx)  # Force in N (unit width)
        end

    end
    
    x_separation = range(start_pos, end_pos, length=n_elem_layered) .* 1e3
    separation = zeros(n_elem_layered)
    for i in 1:n_elem_layered
        elem = start_elem + i - 1
        #separation[i] = abs(u_last[n_dof_steel + (i-1)*2 + 1] - u_last[(elem-1)*2 + 1]) * 1e6
        separation[i] = si_deflection[i] - parylene_deflection[i] 
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

plot_results(u_hist, damage)
    





