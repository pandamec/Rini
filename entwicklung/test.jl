
using Pkg
using Rini
using Plots

####### Fatigue-Test

    Fig=plot()
    Figzoom=plot()

    ## Geometrie 
    E= [150000; 3200; 3200; 200000] #First layer Silicon
    Ge=[0.5 30 15 ; 0.01 30 15; 0.01 30 15; 0.5 15 70 ] 

    sigma1_max=[]
    sigma2_max=[]
    sigma3_max=[]
    sigmas_max=[]

    deltamax=range(0,10,10)

    for deltav in deltamax

        ## Neutral axis berechnung
        sigma = Rini.StaticBeam(E,Ge,deltav)

        push!(sigma3_max,sigma[1])
        push!(sigma2_max,sigma[2])
        push!(sigma1_max,sigma[3])
        push!(sigmas_max,sigma[4])

    end

    plot!(Fig,deltamax, sigma3_max, label="Si", lw=2, linestyle=:dash, color=:black, marker=:square)
    plot!(Fig,deltamax, sigma2_max, label="Parylene 2", lw=2, linestyle=:dash, color=:yellow, marker=:circle)
    plot!(Fig,deltamax, sigma1_max, label="Parylene 1", lw=2, linestyle=:dash, color=:orange, marker=:diamond)
    plot!(Fig,deltamax, sigmas_max, label="Stahl", xlabel="Verschiebungsamplitude(mm)", ylabel="Belastung (MPa)", lw=2, linestyle=:dash, color=:gray, marker=:circle)

    plot!(Figzoom,deltamax, sigma2_max, label="Parylene 2", lw=2, linestyle=:dash, color=:yellow, marker=:circle)
    plot!(Figzoom,deltamax, sigma1_max, label="Parylene 1", xlabel="Verschiebungsamplitude(mm)", ylabel="Belastung (MPa)", lw=2, linestyle=:dash, color=:orange, marker=:diamond)


    Fig
    Figzoom
