
using SymPy
using Plots


#### Strain rate ####

Fig=plot()
Figzoom=plot()

FigStrain=plot()


include("D:/01 Projekt/03 Programm/Rini aktuell/Rini/src/Rini.jl")




struct Material
    E::Float64          # Young's modulus (Pa)
    ν::Float64          #   Poisson's ratio
end

struct BeamSection
    w::Float64          #   Young's modulus (Pa)
    t::Float64          #   Poisson's ratio
    l::Float64          # Thickness (m)
 
end

const Si        =   Material(130e9, 0.28)
const Parylene  =   Material(2.8e9, 0.4)
const Steel     =   Material(210e9, 0.3)

const Si_section      =   BeamSection(5e-3,675e-6,5e-3)
const Parylene_section      =   BeamSection(2e-3,10e-6,2e-3)
const Steel_section      =   BeamSection(15e-3,200e-6,60e-3)


deltamax=range(0,2e-3,10) #mm



    ## Mechaninische Eigenschaften
    #wp=10 #mm
    # Layer 3 Silicon
    #E3=150000
    #Sigma_B3=170 #MPa

    # Layer 2 Parylene
    #E2=3200
    #Sigma_B2=69 #MPa

    # Layer 1 Parylene
    #E1=3200
    #Sigma_B1=69 #MPa

    sigmaSi=[]
    sigmaPa_max=[]
    sigmaPa_min=[]
    sigmaSteel=[]

    c=symbols("c")
    eq=(1/2)*(Steel_section.w*Steel.E*((Steel_section.t-c)^2-(c)^2)+Parylene_section.w*Parylene.E*((Steel_section.t-c+Parylene_section.t)^2-(Steel_section.t-c)^2)+Si_section.w*Si.E*((Steel_section.t-c+Parylene_section.t+Si_section.t)^2-(Steel_section.t-c+Parylene_section.t)^2))

    sol=solve(eq,c)
    sol
    c=sol[]
        


    for i in deltamax

        ## Neutral axis berechnung
        I=(Steel_section.w*Steel_section.t^3)/12
        ## Krafte und Torque
        F=i*48*Steel.E*I/Steel_section.l^3
        M_max=(F/2)*Steel_section.l/2
        x_p=Steel_section.l/2
        M=(x_p/Steel_section.l)*M_max #N.mm



        ## Belastung

        m=symbols("m")
        eqb=M-m*(1/3)*(Steel_section.w*Steel.E*((Steel_section.t-c)^3-(c)^3)+Parylene_section.w*Parylene.E*((Steel_section.t-c+Parylene_section.t)^3-(Steel_section.t-c)^3)+Si_section.w*Si.E*((Steel_section.t-c+Parylene_section.t+Si_section.t)^3-(Steel_section.t-c+Parylene_section.t)^3)+Steel_section.w*Steel.E*((c)^3))
        solb=solve(eqb,m)
        solb
        m=solb[]

        ymax=(Steel_section.t-c)+Parylene_section.t+Si_section.t
        emax=ymax*m
        sigma_Si=m*(Steel_section.t-c+Parylene_section.t+Si_section.t)*Si.E
        sigma_Pa_max=m*(Steel_section.t-c+Parylene_section.t)*Parylene.E
        sigma_Pa_min=m*(Steel_section.t-c)*Steel.E
        sigma_Steel=m*(c)*Steel.E


        push!(sigmaSi,sigma_Si)
        push!(sigmaPa_max,sigma_Pa_max)
        push!(sigmaPa_min,sigma_Pa_min)
        push!(sigmaSteel,sigma_Steel)

    end

    plot!(Fig,deltamax, sigmaSi, label="Si", lw=2, linestyle=:dash, color=:black, marker=:square)
    plot!(Fig,deltamax, sigmaPa_max, label="Parylene max", lw=2, linestyle=:dash, color=:yellow, marker=:circle)
    plot!(Fig,deltamax, sigmaPa_min, label="Parylene min", xlabel="Verschiebungsamplitude(mm)", ylabel="Belastung (MPa)", lw=2, linestyle=:dash, color=:gray, marker=:circle)
    plot!(Fig,deltamax, sigmaSteel, label="Stahl min", xlabel="Verschiebungsamplitude(mm)", ylabel="Belastung (MPa)", lw=2, linestyle=:dash, color=:gray, marker=:circle)

    plot!(FigStrain,deltamax, sigmaSi/Si.E, label="Si", lw=2, linestyle=:dash, color=:black, marker=:square)
    plot!(FigStrain,deltamax, sigmaPa_max/Parylene.E, label="Parylene max", lw=2, linestyle=:dash, color=:yellow, marker=:circle)
    plot!(FigStrain,deltamax, sigmaPa_min/Steel.E, label="Parylene min", xlabel="Verschiebungsamplitude(mm)", ylabel="Strain ", lw=2, linestyle=:dash, color=:gray, marker=:circle)
    plot!(FigStrain,deltamax, sigmaSteel/Steel.E, label="Stahl", xlabel="Verschiebungsamplitude(mm)", ylabel="Strain ", lw=2, linestyle=:dash, color=:gray, marker=:circle)

    FigStrainRate=plot()
    FigSpeedRate=plot()

    freq=1.2

    plot!(FigStrainRate,deltamax, sigmaSi*(4*freq)/Si.E, label="Si", lw=2, linestyle=:dash, color=:black, marker=:square)
    plot!(FigStrainRate,deltamax, sigmaPa_max*(4*freq)/Parylene.E, label="Parylene max", lw=2, linestyle=:dash, color=:yellow, marker=:circle)
    plot!(FigStrainRate,deltamax, sigmaPa_min*(4*freq)/Steel.E, label="Parylene min", xlabel="Verschiebungsamplitude(mm)", ylabel="Strain Rate", lw=2, linestyle=:dash, color=:gray, marker=:circle)
    plot!(FigStrainRate,deltamax, sigmaSteel*(4*freq)/Steel.E, label="Stahl", xlabel="Verschiebungsamplitude(mm)", ylabel="Strain Rate", lw=2, linestyle=:dash, color=:gray, marker=:circle)

    plot!(FigSpeedRate,deltamax*1000, -sigmaPa_min*(4*freq)*60*10/Parylene.E, label="Parylene", xlabel="Verschiebungsamplitude(mm)", ylabel="Speed Rate @10mm, f1.2 (mm/min)",  lw=2, linestyle=:dash, color=:yellow, marker=:circle)
    
    

    Fig
    Figzoom
    FigStrain
    FigSpeedRate