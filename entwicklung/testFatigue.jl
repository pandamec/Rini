
# Extra layer added Klebstoff
#using Pkg
using Rini
using Plots

####### Fatigue-Test

    Fig=plot()
    Figzoom=plot()

    ## Geometrie 
    E= [150000; 2800; 2800; 3.3; 200000] #First layer Silicon
    Ge=[0.32 20 10 ; 0.01 20 10; 0.01 20 10; 0.5 20 10; 0.5 65 15 ] 

    sigma1_max=[]
    sigma2_max=[]
    sigma3_max=[]
    sigmak_max=[]
    sigmas_max=[]

    deltamax=range(0,1,10)

    for deltav in deltamax

        ## Neutral axis berechnung
        sigma = Rini.StaticBeam(E,Ge,deltav)

        push!(sigma3_max,sigma[1])
        push!(sigma2_max,sigma[2])
        push!(sigma1_max,sigma[3])
        push!(sigmak_max,sigma[4])
        push!(sigmas_max,sigma[5])

    end
    Es=200000
    w=Ge[5,3]
    ds=Ge[5,1]
    l=Ge[5,2]
    I=w*ds^3/12
    k=48*Es*I/l^3   

    plot!(Fig,deltamax*k, sigma3_max, label="Si", lw=2, linestyle=:dash, color=:black, marker=:square)
    plot!(Fig,deltamax*k, sigma2_max, label="Parylene 2", lw=2, linestyle=:dash, color=:yellow, marker=:circle)
    plot!(Fig,deltamax*k, sigma1_max, label="Parylene 1", lw=2, linestyle=:dash, color=:orange, marker=:diamond)
    plot!(Fig,deltamax*k, sigmak_max, label="Klebstoff", lw=2, linestyle=:dash, color=:red, marker=:square)
    plot!(Fig,deltamax*k, sigmas_max, label="Stahl", xlabel="Verschiebungsamplitude(mm)", ylabel="Vorspannung (MPa)", lw=2, linestyle=:dash, color=:gray, marker=:circle)

    
    plot!(Figzoom,deltamax*k, sigma2_max, label="Parylene 2", lw=2, linestyle=:dash, color=:yellow, marker=:circle)
    plot!(Figzoom,deltamax*k, sigma1_max, label="Parylene 1",lw=2, linestyle=:dash, color=:orange, marker=:diamond)
    plot!(Figzoom,deltamax*k, sigmak_max, label="Klebstoff", xlabel="Verschiebungsamplitude(mm)", ylabel="Vorspannung (MPa)", lw=2, linestyle=:dash, color=:red, marker=:square)


    Fig
    Figzoom

#### Thermal expansion #####

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
    plot!(FigThermal3, deltaTRange, (dvm/vm0)*100,  ylims=(0, 2),   label="Aluminum"     ,  xlabel="Temperature (C)", ylabel="Volume increase(%)", lw=2, linestyle=:dash, color=:gray, marker=:circle)
    plot!(FigThermal3, deltaTRange, (dvs/vs0)*100,   ylims=(0, 2),   label="Parylene"     ,  xlabel="Temperature (C)", ylabel="Volume increase(%)",lw=2, linestyle=:dash, color=:orange, marker=:square)
