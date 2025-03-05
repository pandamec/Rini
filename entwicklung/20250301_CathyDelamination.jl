
## 05.03.2025 Ver 0.1 Modell-vorschlag
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
    G_c:: Float64   # Critical energy release rate (J/m2)

    function CohesiveProperties(σ_max::Float64, δ_c::Float64, m::Float64)
        K0  = σ_max / δ_c
        G_c = σ_max* δ_c/2
        return new(K0, σ_max, δ_c, m,G_c)
    end
end

struct TestSetup
    L_steel::Float64  
    L_layered::Float64 
    width_layered::Float64
    width_steel::Float64  
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
const Si        =   Material(160e9, 0.28, 675e-6)
const Parylene  =   Material(2.8e9, 0.4, 10e-6)
const Steel     =   Material(210e9, 0.3, 200e-6)


# Geometrie der Probe
const L_steel   = 60e-3
const L_layered = 5e-3
const width_layered=2e-3
const width_steel=15e-3
const n_elem    = 100               # Finite Elemente Methode benutzt fuer die Modellierung
const dx        = L_steel / n_elem
const nodes     = collect(0:dx:L_steel)

const start_pos         = 20e-3
const end_pos           = start_pos + L_layered
const start_elem        = floor(Int, start_pos / dx) + 1
const end_elem          = floor(Int, end_pos / dx)
const n_elem_layered    = end_elem - start_elem + 1
const setup             = TestSetup(L_steel,L_layered,width_layered,width_steel,n_elem,dx,nodes,start_pos,end_pos,start_elem,end_elem,n_elem_layered)

    
#println("Starting simulation...")
#u_hist, damage = simulate_fatigue(setup,max_force, cycles,Si,Parylene,Steel,CZM)
#println("Simulation completed. Generating plots with Makie...")


function update_plot_results!(u_hist, damage, m,damage_history)

    forces = []
    n_dof_steel = length(nodes) * 2
    u_last = u_hist[end]
    u_mid=u_hist[17000]
    x_steel = nodes .* 1e3
    steel_deflection = u_last[1:2:n_dof_steel] .* 1e6
    
    x_layered = range(start_pos, end_pos, length=n_elem_layered + 1) .* 1e3
    parylene_deflection = zeros(n_elem_layered + 1)
    si_deflection = zeros(n_elem_layered + 1)

    # Interface Forces
    for i in 1:n_elem_layered + 1
        elem = start_elem + i - 1
        parylene_deflection[i] = u_last[(elem-1)*2 + 1] * 1e6 - steel_deflection[start_elem]
        si_deflection[i] = u_last[n_dof_steel + (i-1)*2 + 1] * 1e6

        if i <= n_elem_layered
            δ = abs(u_mid[n_dof_steel + (i-1)*2 + 1] - u_mid[(elem-1)*2 + 1])
            push!(forces, CZM.K0 * (1 - damage_history[17000][i]) * δ * dx)
        end
    end
    
    x_separation = range(start_pos, end_pos, length=n_elem_layered) .* 1e3
    separation = zeros(n_elem_layered)
    for i in 1:n_elem_layered
        elem = start_elem + i - 1
        separation[i] = si_deflection[i] - parylene_deflection[i] 
    end
    
    s = []
    for u in u_hist
        si_deflection_last = u[(start_elem+n_elem_layered)*2 + 1] * 1e6 - steel_deflection[start_elem]
        parylene_deflection_last = u[n_dof_steel + (n_elem_layered)*2 + 1] * 1e6
        push!(s, abs(si_deflection_last - parylene_deflection_last))
    end
    crack_over_cycles = s
    cycles = 1:length(u_hist)
    
    # Print statements remain the same
    println("Max separation (μm): ", maximum(separation))
    println("Max damage: ", maximum(damage))
    println("Max steel deflection (μm): ", maximum(abs.(steel_deflection)))
    # ... rest of print statements ...

    # Update plots with new data - ensure labels are always included
    lines!(ax1, x_steel, steel_deflection, color=:green, linewidth=2, label="Steel")
    #lines!(ax1, [start_pos * 1e3, end_pos * 1e3], [0, 0], color=:black, linestyle=:dash, label="Layered Region")  # Changed vlines to lines for legend compatibility

    lines!(ax2, x_layered, parylene_deflection, color=:orange, linewidth=2, label="Parylene: $m")
    lines!(ax2, x_layered, si_deflection, color=:purple, linewidth=2, label="Silicon: $m")
    lines!(ax3, x_separation, damage, linewidth=2, label=m)
    lines!(ax4, cycles, crack_over_cycles, linewidth=2, label=m)
    lines!(ax5, x_separation, forces, linewidth=2, label=m)

    # Update layout
    fig[1:3, 1:2] = [ax1 ax3; ax2 ax4; ax5 GridLayout()]
    return Dict(
        "steel_deflection" => steel_deflection,
        "parylene_deflection" => parylene_deflection,
        "si_deflection" => si_deflection,
        "forces" => forces,
        "separation" => separation,
        "crack_over_cycles" => crack_over_cycles
    )
end


# Variables study 
#m=[ 2 4 6 8 10]*1e-6
#label=[ "m=2e-6" "m=4e-6" "m=6e-6" "m=8e-6" "m=10e-6" ] 

#m=[ 3 5 10 20]*1e-10
m=[  1 ]*1e-8
#label=[  "a" "b" "c" "d" ] 
label=[  "m=6e-6"  ] 

fig = Figure(resolution=(1200, 1000))
ax1 = Axis(fig[1, 1], title="Steel Beam Deformation (3-Point Bending)", xlabel="Position (mm)", ylabel="Deflection (μm)")
vlines!(ax1, [start_pos * 1e3, end_pos * 1e3], color=:black, linestyle=:dash, label="Layered Region")
ax2 = Axis(fig[2, 1], title="Two-Layered Beam Deformation (20–25 mm)", xlabel="Position (mm)", ylabel="Deflection (μm)")
ax3 = Axis(fig[1, 2], title="Damage Distribution", xlabel="Position (mm)", ylabel="Damage")
ax4 = Axis(fig[2, 2], title="Separation at the interface", xlabel="Cycle", ylabel="Separation (μm)")
ax5 = Axis(fig[3, 1:2], title="Internal Forces Between Si and Parylene at 15000 Cycles", xlabel="Position (mm)", ylabel="Force (N)")

# Test setup
max_force = -0.1
cycles = 20000

CZM             =   CohesiveProperties(0.5e6, 0.0025,6e-6)

for i in 1:length(m)
    CZM_i      =   CohesiveProperties(0.5e6, 0.00025, m[i])

    u_hist, damage,damage_history = simulate_fatigue(setup,max_force, cycles,Si,Parylene,Steel,CZM_i)
    
    update_plot_results!(u_hist, damage, label[i],damage_history)
    print("läuft noch")
end

axislegend(ax3)
axislegend(ax4)
axislegend(ax5)