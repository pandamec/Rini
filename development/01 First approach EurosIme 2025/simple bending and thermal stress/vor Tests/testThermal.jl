

using Rini
using Plots

#### Parylene - Aluminium #####

    ###### Cylindrical via ############

    ## Materials
    E=[2800 0.4 35*10^(-6); 70000 0.33 23.6*10^(-6)] #Parylene und Aluminium

    ## Geometrie
    Ge=[0.035 0.030 10] #Substate Metal Lange

    ## Load

    deltaTRange=range(0, 200, step=10) 
    sigma=[]
    dr=[]
    dvm=[]
    dvs=[]

    r_s2, r_m, l  = Ge

    vm0  = (pi*l)*(r_m^2)
    vs0  = (pi*l)*(r_s2^2-r_m^2)

    for deltaT in deltaTRange

        ## Kompatibilität
        E[2,1]=Rini.E_Aluminium(deltaT)

        p, u_m, u_s2= Rini.MetalSubstrat(E,Ge,deltaT)

        dvm_i = (pi*l)*((r_m+u_m)^2-r_m^2)
        dvs_i = (pi*l)*(((r_s2+u_s2)^2-(r_m+u_m)^2)-(r_s2^2-r_m^2))

        push!(sigma, p)
        push!(dr,  u_m)
        push!(dvm, dvm_i)
        push!(dvs, dvs_i)

    end

    FigThermal1=plot()
    FigThermal2=plot()
    FigThermal3=plot()

    plot!(FigThermal1, deltaTRange, dr*1000,    label="u (um)" , lw=2, linestyle=:dash, color=:blue, marker=:square)
    plot!(FigThermal2, deltaTRange, sigma, ylims=(0, 1), label="Parylene-Aluminum"       ,  xlabel="Temperature (C)", ylabel="Vorspannung(MPa)", lw=2, linestyle=:dash, color=:black, marker=:square)
    plot!(FigThermal3, deltaTRange, (dvm/vm0)*100,  ylims=(0, 2.5),   label="Aluminum"     ,  xlabel="Temperature (C)", ylabel="Volume increase(%)", lw=2, linestyle=:dash, color=:gray, marker=:circle)
    plot!(FigThermal3, deltaTRange, (dvs/vs0)*100,   ylims=(0, 2.5),   label="Parylene"     ,  xlabel="Temperature (C)", ylabel="Volume increase(%)",lw=2, linestyle=:dash, color=:orange, marker=:square)


#### Polyimide - Copper #####
    
    ## Materials
    E=[3200 0.34 55*10^(-6); 117000 0.3 16.5*10^(-6)] #Parylene und Aluminium

    ## Geometrie
    Ge=[0.035 0.030 10] #Substate Metal Lange

    ## Load

    deltaTRange=range(0, 200, step=10) 
    sigma=[]
    dr=[]
    dvm=[]
    dvs=[]

    r_s2, r_m, l  = Ge

    vm0  = (pi*l)*(r_m^2)
    vs0  = (pi*l)*(r_s2^2-r_m^2)

    for deltaT in deltaTRange

        ## Kompatibilität
        E[2,1]=Rini.E_Copper(deltaT)

        p, u_m, u_s2= Rini.MetalSubstrat(E,Ge,deltaT)

        dvm_i = (pi*l)*((r_m+u_m)^2-r_m^2)
        dvs_i = (pi*l)*(((r_s2+u_s2)^2-(r_m+u_m)^2)-(r_s2^2-r_m^2))

        push!(sigma, p)
        push!(dr,  u_m)
        push!(dvm, dvm_i)
        push!(dvs, dvs_i)

    end

    Fig2Thermal1=plot()
    Fig2Thermal2=plot()
    Fig2Thermal3=plot()

    plot!(Fig2Thermal1, deltaTRange, dr*1000,    label="u (um)" , lw=2, linestyle=:dash, color=:blue, marker=:square)
    plot!(Fig2Thermal2, deltaTRange, sigma, ylims=(0, 1), label="Polyimide-Copper"       ,  xlabel="Temperature (C)", ylabel="Vorspannung(MPa)", lw=2, linestyle=:dash, color=:black, marker=:square)
    plot!(Fig2Thermal3, deltaTRange, (dvm/vm0)*100,  ylims=(0, 2.5),   label="Copper"     ,  xlabel="Temperature (C)", ylabel="Volume increase(%)", lw=2, linestyle=:dash, color=:gray, marker=:circle)
    plot!(Fig2Thermal3, deltaTRange, (dvs/vs0)*100,   ylims=(0, 2.5),   label="Polyimide"     ,  xlabel="Temperature (C)", ylabel="Volume increase(%)",lw=2, linestyle=:dash, color=:orange, marker=:square)


    