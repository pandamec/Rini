
#### Log
    ## 05.03.2025 VAer 0.1 Modell-vorschlag
    ## 07.03.2025 Ver  0.2 Crack growth idea added
    ## 08.03.2025 Ver  0.3 Optimization algorithm started
    ## 09.03.2025 Ver  0.4 Optimisierung-vorschlag. FCZM Funktion getrennt

###############################
using Rini
using GLMakie
# Eigenschaften der Materialien

# Materialien der Layers
const Si        =   Material(130e9, 0.28, 675e-6)
const Parylene  =   Material(2.8e9, 0.4, 10e-6)
const Steel     =   Material(210e9, 0.3, 200e-6)

# Geometrie der Probe
const L_steel   = 60e-3
const L_layered = 5e-3
const width_layered=2e-3
const width_steel=15e-3
const n_elem    = 100              # Finite Elemente Methode benutzt fuer die Modellierung
const dx        = L_steel / n_elem
const nodes     = collect(0:dx:L_steel)

const start_pos         = 20e-3
const end_pos           = start_pos + L_layered
const start_elem        = floor(Int, start_pos / dx) + 1
const end_elem          = floor(Int, end_pos / dx)
const n_elem_layered    = end_elem - start_elem + 1
const setup             = TestSetup(L_steel,L_layered,width_layered,width_steel,n_elem,dx,nodes,start_pos,end_pos,start_elem,end_elem,n_elem_layered)

# Plots
function update_plot_results!(CZM,u_hist, damage, m,damage_history,Gc_history,a_history)

    forces = []
    interface_stress=[]
    n_dof_steel = length(nodes) * 2
    u_last = u_hist[end]
    cycle10=Int(length(u_hist)/10)
    u_mid=u_hist[cycle10]
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
            push!(forces, CZM.K0 * (1 - damage_history[cycle10][i]) * δ * dx)
            push!(interface_stress, CZM.K0 * (1 - damage_history[cycle10][i]))
        
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
    lines!(ax5, x_separation, interface_stress*1e-6, linewidth=2, label=m)
    lines!(ax6, cycles, Gc_history, linewidth=2, label=m)
    lines!(ax7, cycles, a_history,linewidth=2, label=m)
end

function export_plot_results!(CZM,u_hist, damage, m,damage_history,Gc_history,a_history)
    forces = []
    interface_stress=[]
    n_dof_steel = length(nodes) * 2
    u_last = u_hist[end]
    cycle10=Int(length(u_hist)/10)
    u_mid=u_hist[cycle10]
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
            push!(forces, CZM.K0 * (1 - damage_history[cycle10][i]) * δ * dx)
            push!(interface_stress, CZM.K0 * (1 - damage_history[cycle10][i]))
        
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
    
    # Update plots with new data - ensure labels are always included
    lines!(axx1, cycles, Gc_history, linewidth=4, color=:orange, linestyle=:dot)
    exp_points=[(1000,Gc_history[1000]),(2000,Gc_history[2000]),(3000,Gc_history[3000]),(4000,Gc_history[4000]),(5000,Gc_history[5000])]


    for (x_exp, y_exp) in exp_points
        vlines!(axx1, [x_exp],color=:black, linestyle=:dash, linewidth=1)
        hlines!(axx1, [y_exp], color=:black, linestyle=:dash, linewidth=1)
    end

end

function exportSeparation_plot_results!(CZM,u_hist, damage, m,damage_history,Gc_history,a_history)
    forces = []
    interface_stress=[]
    n_dof_steel = length(nodes) * 2
    u_last = u_hist[end]
    cycle10=Int(length(u_hist)/10)
    u_mid=u_hist[cycle10]
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
            push!(forces, CZM.K0 * (1 - damage_history[cycle10][i]) * δ * dx)
            push!(interface_stress, CZM.K0 * (1 - damage_history[cycle10][i]))
        
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
    
    # Update plots with new data - ensure labels are always included
    lines!(axx2, cycles, crack_over_cycles, linewidth=4, color=:blue, linestyle=:dot)
    exp_points=[(1000,crack_over_cycles[1000]),(2000,crack_over_cycles[2000]),(3000,crack_over_cycles[3000]),(4000,crack_over_cycles[4000]),(5000,crack_over_cycles[5000])]


    for (x_exp, y_exp) in exp_points
        vlines!(axx2, [x_exp],color=:black, linestyle=:dash, linewidth=1)
        hlines!(axx2, [y_exp], color=:black, linestyle=:dash, linewidth=1)
    end

end


########### Optimization algorithm ###############

function main(CZM0,fatigueExp,max_force,N_cycles)
    
    # Run optimization
    optimized_params, error = Rini.optimize_FCZM(CZM0,setup,max_force,N_cycles,fatigueExp; max_iterations=n_it,Si,Parylene,Steel)
    
    # Print results
    println("Optimized Parameters:")
    println("σ_max = ", optimized_params[1])
    println("δ_max = ", optimized_params[2])
    println("δ_c, = ", optimized_params[3])
    println("m, = ", optimized_params[4])
    println("error = ", error)
    CZM_fit = CohesiveProperties(optimized_params[1],optimized_params[2],optimized_params[3],optimized_params[4])
    
    return CZM_fit
end

############## Example 

# Test Bedingungen
max_force       = -0.1
N_cycles        =  5000  # Number of cycles for the experiment
#CZM0            =   CohesiveProperties(0.3e5, 0.0025/3,0.0025,1e-8) #Initial parameters of the Cohesive Zone Model
fatigueExp      =   FatigueData(44,28e-6) # Example experimental data 
n_it            = 1

#CZM0=CohesiveProperties( 0.1e6, 0.0004, 0.002,1e-8)
CZM0=CohesiveProperties(1e6, 0.00025/2,0.0002,0.2e-7)

CZM_fit=main(CZM0,fatigueExp,max_force,N_cycles) # Find the better adjusment for the CZM parameter
CZM_fit= CohesiveProperties(999999.9999999994, 0.0001150001830588302, 0.00017444895389732032, 1.9999998230835003e-8)

############### Plots Simulation Fit ###########

label="model"

    fig = Figure(resolution=(1200, 1000))

    #ax1 = Axis(fig[1, 1], title="Steel Beam Deformation (3-Point Bending)", xlabel="Position (mm)", ylabel="Deflection (μm)")
    #vlines!(ax1, [start_pos * 1e3, end_pos * 1e3], color=:black, linestyle=:dash, label="Layered Region")
    #ax2 = Axis(fig[2, 1], title="Two-Layered Beam Deformation (20–25 mm)", xlabel="Position (mm)", ylabel="Deflection (μm)")
    #ax3 = Axis(fig[1, 2], title="Damage Distribution", xlabel="Position (mm)", ylabel="Damage")
    #ax4 = Axis(fig[2, 2], title="Separation at the interface", xlabel="Cycle", ylabel="Separation (μm)")
    #ax5 = Axis(fig[3, 1], title="Internal Stress Between Si and Parylene at 10% Cycles", xlabel="Position (mm)", ylabel="Stress (MPa)")
    #ax6  = Axis(fig[3, 2], title="Critical energy release rate Gc over cycles", xlabel="Cycle", ylabel="Gmax (J/m2)")
    #ax7  = Axis(fig[4, 1], title="Crack growth", xlabel="Cycle", ylabel="a (um)")

    ax1     = Axis(fig[1, 1], xlabel="Position (mm)", ylabel="Deflection (μm)")
    vlines!(ax1, [start_pos * 1e3, end_pos * 1e3], color=:black, linestyle=:dash, label="Layered Region")
    ax2     = Axis(fig[2, 1], xlabel="Position (mm)", ylabel="Deflection (μm)",xgridvisible = false,ygridvisible = false)
    ax3     = Axis(fig[1, 2], xlabel="Position (mm)", ylabel="Damage",xgridvisible = false,ygridvisible = false)
    ax4     = Axis(fig[2, 2], xlabel=L"Cycle", ylabel=L"\delta \, (\mu m)", xgridvisible = false, ygridvisible = false)
    ax5     = Axis(fig[3, 1], xlabel=L"Position (mm)", ylabel="Stress (MPa)",xgridvisible = false,ygridvisible = false)
    ax6     = Axis(fig[3, 2], xlabel=L"Cycle", ylabel=L"G_{\text{max}} \, (\text{J/m}^2)",xgridvisible = false,ygridvisible = false)
    ax7     = Axis(fig[4, 1], xlabel=L"Cycle", ylabel=L"a \, (\mu m)",xgridvisible = false,ygridvisible = false,xlabelcolor = :black,
    ylabelcolor = :black)


    u_hist, damage,damage_history,Gc_history,a_history, CZM_hist = simulate_fatigue(setup,max_force, N_cycles,Si,Parylene,Steel,CZM_fit)
            
    update_plot_results!(CZM_fit,u_hist, damage, label,damage_history,Gc_history,a_history)
    println("Gc (J/m2): ", CZM_fit.G_c)

    axislegend(ax3)
    axislegend(ax4)
    axislegend(ax5)
    axislegend(ax6)
    axislegend(ax7)
    
    fig2 = Figure()  # High resolution
    axx1 = Axis(fig2[1, 1],
     xlabel=L"Cycle", 
     ylabel=L"G_{\text{max}} \, (\text{J/m}^2)",
     xgridvisible = false,ygridvisible = false,
     xticks=0:1000:5000,
     yticks=60:2:90,
     xlabelsize=24,
      ylabelsize=24,
      xticklabelsize=20,
       yticklabelsize=20)
    export_plot_results!(CZM_fit,u_hist, damage, label,damage_history,Gc_history,a_history)
    
    fig3 = Figure()  # High resolution
    axx2 = Axis(fig3[1, 1],
     xlabel=L"Cycle", 
     ylabel=L"\delta \, (\mu m)",
     xgridvisible = false,ygridvisible = false,
     xticks=0:1000:5000,
     yticks=0:0.1:25,
     xlabelsize=24,
      ylabelsize=24,
      xticklabelsize=20,
       yticklabelsize=20)
    exportSeparation_plot_results!(CZM_fit,u_hist, damage, label,damage_history,Gc_history,a_history)

###################### Cycle evaluation ########

dN=5000

CZM_i=Rini.FCZM(CZM_fit,0,dN,0)
CZM_fit
