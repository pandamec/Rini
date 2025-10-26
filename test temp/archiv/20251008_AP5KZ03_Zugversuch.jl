using Pkg
using CSV
using DataFrames
using Polynomials
using CairoMakie

## Geometrie
w=10*1e-3
t=15*1e-6
l=75*1e-3

As=w*t

using Plots

# Plot the first curve with lines and markers
fig = Figure(resolution = (1200, 1200))
font=24
ax1 = Axis(fig[1,1],
    xlabel = L"\epsilon (%)",
    ylabel = L"\sigma (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false
)


ax2 = Axis(fig[1,2],
    xlabel = L"\epsilon (%)",
    ylabel = L"\sigma (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false
)

ax3 = Axis(fig[1,3],
    xlabel = L"\epsilon (%)",
    ylabel = L"\sigma (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false
)


ax4 = Axis(fig[2,1],
    xlabel = L"v (%/min)",
    ylabel = L"E (GPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false
)

ax5 = Axis(fig[2,2],
    xlabel = L"v (%/min)",
    ylabel = L"\sigma_B (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false
)

ax6 = Axis(fig[2,3],
    xlabel = L"v (%/min)",
    ylabel = L"\sigma_F (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false
)


function import_TIRA(file_path)
    
    df=CSV.read(file_path,DataFrame) 
    df.Zeit = parse.(Float64, replace.(df.Zeit, "," => "."))
    df.Länge = parse.(Float64, replace.(df.Länge, "," => "."))
    df.Weg = parse.(Float64, replace.(df.Weg, "," => "."))
    #df.Dehnung = parse.(Float64, replace.(df.Dehnung, "," => "."))
    df.Kraft = parse.(Float64, replace.(df.Kraft, "," => "."))
    df.dL_ORG = parse.(Float64, replace.(df.dL_ORG, "," => "."))

    return df

end

function filter_range(df::DataFrame, column::String, low, high)
    return filter(row -> low ≤ row[column] ≤ high, df)
end


#### Probe 1

    # Replace with the ath to your Excel file


    file_path="test/AP5ZK03/251006-01.csv"
    df=import_TIRA(file_path)

    df_filtered=filter_range(df,"dL_ORG",0,5)
    T=df_filtered[!,:Zeit]
    X=df_filtered[!,:dL_ORG]
    Y=df_filtered[!,:Kraft]


    # Fit a 1st-degree polynomial (linear)
    p_strain=fit(T[1:1000],X[1:1000]*60,1)

    # Fit a 1st-degree polynomial (linear)
    p_strain=fit(X[500:2000],(Y[500:2000]*1e-6)/As,1)
    Xfit=X[500:2000]
    Yfit = p_strain.(Xfit)

    # Fit a 1st-degree polynomial (linear)
    p_strain=fit(X[500:2000],(Y[500:2000]*1e-6)/As,1)
    Xfit=X[500:2000]
    Yfit = p_strain.(Xfit)


    # Plot the first curve with lines and markers


    Makie.scatter!(ax1,X,Y*1e-6/As; markersize=12, label="Strain rate 1.25 %/min",color=:orange)
    #Makie.lines!(ax,Xfit, Yfit; linewidth=5, label= "linear fit",color=:black)
    axislegend(ax1, position=:rb, labelsize=font-10, framevisible=false)

    #save("tensile.svg",fig)


#### Gesamtgruppe Raum T

    BaseName="test/AP5ZK03/251006-0"
    name="251006-0"
    df_Group=[]

    for i in [7,8,1,2,3,4,5,6]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0,3)
        push!(df_Group,df_filtered)
    end


    # Plot the first curve with lines and markers

    palette = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :black]


    log=[]
    log_total=[]
    for i in df_Group
        
        T=i[!,:Zeit]
        X=i[!,:dL_ORG]
        Y=i[!,:Stress]
        df_linear=filter_range(i,"Stress",0,20)
        df_linear=filter_range(df_linear,"dL_ORG",0.1,0.25)
        T_linear=df_linear[!,:Zeit]
        X_linear=df_linear[!,:dL_ORG]
        Y_linear=df_linear[!,:Stress]

        sigma0_25=maximum(Y_linear)
        sigmaB=maximum(Y)
        sigmaF=maximum(Y_linear)
        # Fit a 1st-degree polynomial (linear)
        fit_X_T=fit(T_linear,X_linear*60,1)
        
        strainRate=fit_X_T[1]
        
        # Fit a 1st-degree polynomial (linear)
        p_fit=fit(X_linear/100,Y_linear,1)
        E_modul=p_fit[1]/1000

        Makie.scatter!(ax2,X,Y; markersize=12, label = "$(round(strainRate, digits=3)) [%/min]")
        
        Makie.scatter!(ax3,X_linear,Y_linear; markersize=12, label = "$(round(strainRate, digits=3)) [%/min]")
        
        #Makie.lines!(ax,Xfit, Yfit; linewidth=5, label= "linear fit",color=:black)
        #axislegend(ax2, labelsize=font-10, framevisible=false)
        #axislegend(ax3, labelsize=font-10, framevisible=false)
        sum=Dict("Probe"=>i[1,:Name],
            " Strain Rate[%/min]"=>strainRate,
                "E[GPa] "=>E_modul,
                "sigma0_25% "=> sigma0_25,
                "sigmaMax"=> sigmaB,
                "Temperature [C]"=> 22)
        push!(log,sum)
        push!(log_total,sum)
        
    end

    axislegend(ax2, labelsize=font-10, framevisible=false)
    axislegend(ax3, labelsize=font-10, framevisible=false)

    df_log=DataFrame(log)

    epsilon=df_log[!," Strain Rate[%/min]"]
    E=df_log[!,"E[GPa] "]
    Makie.scatter!(ax4,epsilon,E)
    Makie.lines!(ax4,epsilon,E)

    Makie.scatter!(ax5,df_log[!," Strain Rate[%/min]"],df_log[!,"sigmaMax"])
    Makie.lines!(ax5,df_log[!," Strain Rate[%/min]"],df_log[!,"sigmaMax"])


    Makie.scatter!(ax6,df_log[!," Strain Rate[%/min]"],df_log[!,"sigma0_25% "])
    Makie.lines!(ax6,df_log[!," Strain Rate[%/min]"],df_log[!,"sigma0_25% "])
    fig


#### Probe unter 50 °C 


    # Plot the first curve with lines and markers
    fig2 = Figure(resolution = (1200, 1200))
    font=24
    ax1_2 = Axis(fig2[1,1],
        xlabel = L"\epsilon (%)",
        ylabel = L"\sigma (MPa)",
        xlabelsize = font,
        ylabelsize = font,
        xticklabelsize = font - 2,
        yticklabelsize = font - 2,
        xgridstyle = :dash,        # dashed grid
        ygridstyle = :dash,
        xgridvisible = false,
        ygridvisible = false
    )


        ax2_2 = Axis(fig2[1,2],
        xlabel = L"\epsilon (%)",
        ylabel = L"\sigma (MPa)",
        xlabelsize = font,
        ylabelsize = font,
        xticklabelsize = font - 2,
        yticklabelsize = font - 2,
        xgridstyle = :dash,        # dashed grid
        ygridstyle = :dash,
        xgridvisible = false,
        ygridvisible = false
    )

    ax3_2 = Axis(fig2[1,3],
    xlabel = L"\epsilon (%)",
    ylabel = L"\sigma (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false
    )

    ax4_2 = Axis(fig2[1,4],
    xlabel = L"\epsilon (%)",
    ylabel = L"\sigma (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false
    )


    ax5_2 = Axis(fig2[2,1],
        xlabel = L"v (%/min)",
        ylabel = L"E (GPa)",
        xlabelsize = font,
        ylabelsize = font,
        xticklabelsize = font - 2,
        yticklabelsize = font - 2,
        xgridstyle = :dash,        # dashed grid
        ygridstyle = :dash,
        xgridvisible = false,
        ygridvisible = false
    )

    ax6_2 = Axis(fig2[2,2],
        xlabel = L"v (%/min)",
        ylabel = L"\sigma_B (MPa)",
        xlabelsize = font,
        ylabelsize = font,
        xticklabelsize = font - 2,
        yticklabelsize = font - 2,
        xgridstyle = :dash,        # dashed grid
        ygridstyle = :dash,
        xgridvisible = false,
        ygridvisible = false
    )

    ax7_2 = Axis(fig2[2,3],
        xlabel = L"v (%/min)",
        ylabel = L"\sigma_F (MPa)",
        xlabelsize = font,
        ylabelsize = font,
        xticklabelsize = font - 2,
        yticklabelsize = font - 2,
        xgridstyle = :dash,        # dashed grid
        ygridstyle = :dash,
        xgridvisible = false,
        ygridvisible = false
    )

    
    file_path="test/AP5ZK03/251010-01.csv"
    df=import_TIRA(file_path)    
    df[!,:Stress]=df[!,:Kraft]*1e-6/As
    df[!,:Strain]=df[!,:dL_ORG]/100
    
    Makie.scatter!(ax1_2,df[!,:Strain],df[!,:Stress]; markersize=12, label="20.6 %/min",color=:orange)
    #Makie.lines!(ax,Xfit, Yfit; linewidth=5, label= "linear fit",color=:black)
    axislegend(ax1_2, position=:rb, labelsize=font-10, framevisible=false)

    df_filtered=filter_range(df,"dL_ORG",0.1,5)
    T=df_filtered[!,:Zeit]
    X=df_filtered[!,:dL_ORG]
    Y=df_filtered[!,:Stress]

    df_linear=filter_range(df_filtered,"dL_ORG",0.1,0.5)
    T_linear=df_linear[!,:Zeit]
    X_linear=df_linear[!,:dL_ORG]
    Y_linear=df_linear[!,:Stress]


    # Fit a 1st-degree polynomial (linear)
    p_strain=fit(T_linear,X_linear*60,1)

    # Fit a 1st-degree polynomial (linear)
    p_strain=fit(X_linear/100,Y_linear,1)

    # Plot the first curve with lines and markers

    Makie.scatter!(ax2_2,X,Y; markersize=12, label="20.6 %/min",color=:orange)
    #Makie.lines!(ax,Xfit, Yfit; linewidth=5, label= "linear fit",color=:black)
    axislegend(ax2_2, position=:rb, labelsize=font-10, framevisible=false)

    fig2
    #save("tensile.svg",fig)

#### Gesamtgruppe um 50 °C

    BaseName="test/AP5ZK03/251010-0"
    name="251010-0"
    df_Group=[]

   for i in [1,2]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0.1,5)
        push!(df_Group,df_filtered)
    end


    # Plot the first curve with lines and markers

    palette = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :black]


    log=[]
    for i in df_Group
        
        T=i[!,:Zeit]
        X=i[!,:dL_ORG]
        Y=i[!,:Stress]
        #df_linear=filter_range(i,"Stress",0,20)
        df_linear=filter_range(i,"dL_ORG",0.3,1)
        T_linear=df_linear[!,:Zeit]
        X_linear=df_linear[!,:dL_ORG]
        Y_linear=df_linear[!,:Stress]

        sigma0_25=maximum(Y_linear)
        sigmaB=maximum(Y)
        sigmaF=maximum(Y_linear)
        # Fit a 1st-degree polynomial (linear)
        fit_X_T=fit(T_linear,X_linear*60,1)
        
        strainRate=fit_X_T[1]
        
        # Fit a 1st-degree polynomial (linear)
        p_fit=fit(X_linear/100,Y_linear,1)
        E_modul=p_fit[1]/1000

        Makie.scatter!(ax3_2,X,Y; markersize=12, label = "$(round(strainRate, digits=3)) [%/min]")
        
        Makie.scatter!(ax4_2,X_linear,Y_linear; markersize=12, label = "$(round(strainRate, digits=3)) [%/min]")
        
        #Makie.lines!(ax,Xfit, Yfit; linewidth=5, label= "linear fit",color=:black)
        #axislegend(ax2, labelsize=font-10, framevisible=false)
        #axislegend(ax3, labelsize=font-10, framevisible=false)
        sum=Dict("Probe"=>i[1,:Name],
            " Strain Rate[%/min]"=>strainRate,
                "E[GPa] "=>E_modul,
                "sigma0_25% "=> sigma0_25,
                "sigmaMax"=> sigmaB,
                "Temperature [C]"=> 50)
        push!(log,sum)
        push!(log_total,sum)
        
    end

    axislegend(ax3_2, labelsize=font-10, framevisible=false)
    axislegend(ax4_2, labelsize=font-10, framevisible=false)

    df_log=DataFrame(log)

    epsilon=df_log[!," Strain Rate[%/min]"]
    E=df_log[!,"E[GPa] "]
    Makie.scatter!(ax5_2,epsilon,E)
    Makie.lines!(ax5_2,epsilon,E)

    Makie.scatter!(ax6_2,df_log[!," Strain Rate[%/min]"],df_log[!,"sigmaMax"])
    Makie.lines!(ax6_2,df_log[!," Strain Rate[%/min]"],df_log[!,"sigmaMax"])


    Makie.scatter!(ax7_2,df_log[!," Strain Rate[%/min]"],df_log[!,"sigma0_25% "])
    Makie.lines!(ax7_2,df_log[!," Strain Rate[%/min]"],df_log[!,"sigma0_25% "])
    fig2



#### Probe unter 100 °C 


    # Plot the first curve with lines and markers
    fig3 = Figure(resolution = (1200, 1200))
    font=24
    ax1_3 = Axis(fig3[1,1],
        xlabel = L"\epsilon (%)",
        ylabel = L"\sigma (MPa)",
        xlabelsize = font,
        ylabelsize = font,
        xticklabelsize = font - 2,
        yticklabelsize = font - 2,
        xgridstyle = :dash,        # dashed grid
        ygridstyle = :dash,
        xgridvisible = false,
        ygridvisible = false
    )


        ax2_3 = Axis(fig3[1,2],
        xlabel = L"\epsilon (%)",
        ylabel = L"\sigma (MPa)",
        xlabelsize = font,
        ylabelsize = font,
        xticklabelsize = font - 2,
        yticklabelsize = font - 2,
        xgridstyle = :dash,        # dashed grid
        ygridstyle = :dash,
        xgridvisible = false,
        ygridvisible = false
    )

    ax3_3 = Axis(fig3[1,3],
    xlabel = L"\epsilon (%)",
    ylabel = L"\sigma (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false
    )

    ax4_3 = Axis(fig3[1,4],
    xlabel = L"\epsilon (%)",
    ylabel = L"\sigma (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false
    )


    ax5_3 = Axis(fig3[2,1],
        xlabel = L"v (%/min)",
        ylabel = L"E (GPa)",
        xlabelsize = font,
        ylabelsize = font,
        xticklabelsize = font - 2,
        yticklabelsize = font - 2,
        xgridstyle = :dash,        # dashed grid
        ygridstyle = :dash,
        xgridvisible = false,
        ygridvisible = false
    )

    ax6_3 = Axis(fig3[2,2],
        xlabel = L"v (%/min)",
        ylabel = L"\sigma_B (MPa)",
        xlabelsize = font,
        ylabelsize = font,
        xticklabelsize = font - 2,
        yticklabelsize = font - 2,
        xgridstyle = :dash,        # dashed grid
        ygridstyle = :dash,
        xgridvisible = false,
        ygridvisible = false
    )

    ax7_3 = Axis(fig3[2,3],
        xlabel = L"v (%/min)",
        ylabel = L"\sigma_F (MPa)",
        xlabelsize = font,
        ylabelsize = font,
        xticklabelsize = font - 2,
        yticklabelsize = font - 2,
        xgridstyle = :dash,        # dashed grid
        ygridstyle = :dash,
        xgridvisible = false,
        ygridvisible = false
    )

    
    file_path="test/AP5ZK03/251006-9.csv"
    df=import_TIRA(file_path)    
    df[!,:Stress]=df[!,:Kraft]*1e-6/As
    df[!,:Strain]=df[!,:dL_ORG]/100
    
    Makie.scatter!(ax1_3,df[!,:Strain],df[!,:Stress]; markersize=12, label="20.6 %/min",color=:orange)
    #Makie.lines!(ax,Xfit, Yfit; linewidth=5, label= "linear fit",color=:black)
    axislegend(ax1_3, position=:rb, labelsize=font-10, framevisible=false)

    df_filtered=filter_range(df,"dL_ORG",0.1,5)
    T=df_filtered[!,:Zeit]
    X=df_filtered[!,:dL_ORG]
    Y=df_filtered[!,:Stress]

    df_linear=filter_range(df_filtered,"dL_ORG",0.1,0.5)
    T_linear=df_linear[!,:Zeit]
    X_linear=df_linear[!,:dL_ORG]
    Y_linear=df_linear[!,:Stress]


    # Fit a 1st-degree polynomial (linear)
    p_strain=fit(T_linear,X_linear*60,1)

    # Fit a 1st-degree polynomial (linear)
    p_strain=fit(X_linear/100,Y_linear,1)

    # Plot the first curve with lines and markers

    Makie.scatter!(ax2_3,X,Y; markersize=12, label="20.6 %/min",color=:orange)
    #Makie.lines!(ax,Xfit, Yfit; linewidth=5, label= "linear fit",color=:black)
    axislegend(ax2_3, position=:rb, labelsize=font-10, framevisible=false)

    fig3
    #save("tensile.svg",fig)

#### Gesamtgruppe um 100 °C

    BaseName="test/AP5ZK03/251006-"
    name="251006-"
    df_Group=[]

   for i in [9,10]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        s0=df[1,:dL_ORG]
        df[!,:dL_ORG]=(df[!,:dL_ORG].-s0)
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0.1,10)
        push!(df_Group,df_filtered)
    end


    # Plot the first curve with lines and markers

    palette = [:red, :blue, :green, :orange, :purple, :cyan, :magenta, :black]


    
    log=[]
    for i in df_Group
        
        T=i[!,:Zeit]
        X=i[!,:dL_ORG]
        Y=i[!,:Stress]
        #df_linear=filter_range(i,"Stress",0,20)
        df_linear=filter_range(i,"dL_ORG",0.5,1.5)
        T_linear=df_linear[!,:Zeit]
        X_linear=df_linear[!,:dL_ORG]
        Y_linear=df_linear[!,:Stress]

        sigma0_25=maximum(Y_linear)
        sigmaB=maximum(Y)
        sigmaF=maximum(Y_linear)
        # Fit a 1st-degree polynomial (linear)
        fit_X_T=fit(T_linear,X_linear*60,1)
        
        strainRate=fit_X_T[1]
        
        # Fit a 1st-degree polynomial (linear)
        p_fit=fit(X_linear/100,Y_linear,1)
        E_modul=p_fit[1]/1000

        Makie.scatter!(ax3_3,X,Y; markersize=12, label = "$(round(strainRate, digits=3)) [%/min]")
        
        Makie.scatter!(ax4_3,X_linear,Y_linear; markersize=12, label = "$(round(strainRate, digits=3)) [%/min]")
        
        #Makie.lines!(ax,Xfit, Yfit; linewidth=5, label= "linear fit",color=:black)
        #axislegend(ax2, labelsize=font-10, framevisible=false)
        #axislegend(ax3, labelsize=font-10, framevisible=false)

        sum=Dict("Probe"=>i[1,:Name],
                " Strain Rate[%/min]"=>strainRate,
                "E[GPa] "=>E_modul,
                "sigma0_25% "=> sigma0_25,
                "sigmaMax"=> sigmaB,
                "Temperature [C]"=> 100)
        push!(log,sum)
        push!(log_total,sum)
        
    end

    axislegend(ax3_3, labelsize=font-10, framevisible=false)
    axislegend(ax4_3, labelsize=font-10, framevisible=false)

    df_log=DataFrame(log)
    df_log_total=DataFrame(log_total)
    epsilon=df_log[!," Strain Rate[%/min]"]
    E=df_log[!,"E[GPa] "]
    Makie.scatter!(ax5_3,epsilon,E)
    Makie.lines!(ax5_3,epsilon,E)

    Makie.scatter!(ax6_3,df_log[!," Strain Rate[%/min]"],df_log[!,"sigmaMax"])
    Makie.lines!(ax6_3,df_log[!," Strain Rate[%/min]"],df_log[!,"sigmaMax"])


    Makie.scatter!(ax7_3,df_log[!," Strain Rate[%/min]"],df_log[!,"sigma0_25% "])
    Makie.lines!(ax7_3,df_log[!," Strain Rate[%/min]"],df_log[!,"sigma0_25% "])
    fig3
