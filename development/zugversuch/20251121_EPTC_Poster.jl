using Pkg
using CSV
using DataFrames
using Polynomials
using CairoMakie
using Statistics


using SymPy
using LinearAlgebra
using SparseArrays
using Statistics
using Optim
    


# Define the predictive function

function import_TIRA(file_path)
    
    df=CSV.read(file_path,DataFrame) 
    df.Zeit = parse.(Float64, replace.(df.Zeit, "," => "."))
    df.Länge = parse.(Float64, replace.(df.Länge, "," => "."))
    df.Weg = parse.(Float64, replace.(df.Weg, "," => "."))
    df.Kraft = parse.(Float64, replace.(df.Kraft, "," => "."))
    df.dL_ORG = parse.(Float64, replace.(df.dL_ORG, "," => "."))

    return df

end

function filter_range(df::DataFrame, column::String, low, high)
    return filter(row -> low ≤ row[column] ≤ high, df)
end


############## Example 

## Geometrie
w=10*1e-3 # Width
t=15*1e-6 # thickness
l=75*1e-3 # length
As=w*t    # stress section

σ_exp=[]
ϵ_exp=[]
dϵdt=[]
##### Strain Rate 1.43 %/min

BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/AP5ZK03/251006-0"
name="251006-"
df_Group=[]
df_Group_Total=[]

for i in [1]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0.1,1)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end

df_Group
push!(σ_exp,df_Group[1][!,:Stress]*1e6)
push!(ϵ_exp,df_Group[1][!,:Strain])
push!(dϵdt,0.0143/60)


##### Strain Rate 4.44 %/min

BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/AP5ZK03/251006-0"
name="251006-"
df_Group=[]
df_Group_Total=[]

for i in [6]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0.1,1)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end

df_Group
push!(σ_exp,df_Group[1][!,:Stress]*1e6)
push!(ϵ_exp,df_Group[1][!,:Strain])
push!(dϵdt,0.0444/60)


### Strain rate 27%
BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/AP5ZK04/251103-"
name="251103-"
df_Group=[]
df_Group_Total=[]

for i in [1]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0.1,1)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end


push!(σ_exp,df_Group[1][!,:Stress]*1e6)
push!(ϵ_exp,df_Group[1][!,:Strain])
push!(dϵdt,0.2703/60)


# Plot the first curve with lines and markers


fig1 = Figure(resolution = (1000, 600));
fig2 = Figure(resolution = (1000, 600));
fig3 = Figure(resolution = (1000, 600));
font=42
ax1_1 = Axis(fig1[1,1],
    xlabel = L"\epsilon (%)",
    ylabel = L"\sigma (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridstyle = :dash,        # dashed grid
    ygridstyle = :dash,
    xgridvisible = false,
    ygridvisible = false,
    limits = ((0, 2), nothing)
)

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


## Fig2 Temperatur Einfluss

σ_exp=[]
ϵ_exp=[]
dϵdt=[]
## Room 27 %/min

BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/AP5ZK04/251103-"
name="251103-"
df_Group=[]

df_Group_Total=[]

for i in [1]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0,2.1)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end

push!(σ_exp,df_Group[1][!,:Stress]*1e6)
push!(ϵ_exp,df_Group[1][!,:Strain])
push!(dϵdt,0.2703/60)

## 50 27 %/min

BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/AP5ZK04/251103-"
name="251103-"
df_Group=[]
df_Group_Total=[]

for i in [7]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0,2.1)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end

push!(σ_exp,df_Group[1][!,:Stress]*1e6)
push!(ϵ_exp,df_Group[1][!,:Strain])
push!(dϵdt,0.2545/60)


## 100 27 %/min

BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/AP5ZK04/251103-"
name="251103-"
df_Group=[]
df_Group_Total=[]

for i in [10]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:dL_ORG]=df[!,:dL_ORG].-0.6
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0,2.1)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end

push!(σ_exp,df_Group[1][!,:Stress]*1e6)
push!(ϵ_exp,df_Group[1][!,:Strain])
push!(dϵdt,0.2428/60)


## Fig2


    Makie.scatter!(ax1_2 ,ϵ_exp[1]*100,σ_exp[1]/1e6;
     markersize=14, 
     label = "T22°C v=$(round(dϵdt[1]*60*100, digits=3)) %/min",
     color = :black)
    
    Makie.scatter!(ax1_2,ϵ_exp[2]*100,σ_exp[2]/1e6;
      label = "T50°C v=$(round(dϵdt[2]*60*100, digits=3)) %/min",
      markersize=14,
      color = :orange)
    
    Makie.scatter!(ax1_2,ϵ_exp[3]*100,σ_exp[3]/1e6;
    label = "T100°C v=$(round(dϵdt[3]*60*100, digits=3)) %/min",
    markersize=14, 
    color = :red)
    axislegend(ax1_2, labelsize=font-10, framevisible=false,position = :rb)

    fig2


### Polyimide Comparison



fig1 = Figure(resolution = (1000, 600));
font=42
ax1_1 = Axis(fig1[1,1],
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




BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/Dharmadasa2020.csv"

df_Polyimide = CSV.read(BaseName, DataFrame;
    delim=';',
    decimal=','
)

df_Polyimide=filter_range(df_Polyimide,"Strain",0,0.02)

## Fig1


    Makie.scatter!(ax1_1 ,ϵ_exp[1]*100,σ_exp[1]/1e6;
     markersize=14, 
     label = "Parylene v=$(round(dϵdt[1]*60*100, digits=3)) %/min ",
     color = :black,
     clip = true)
    
    Makie.lines!(ax1_1,df_Polyimide[!,1]*100,df_Polyimide[!,2].*140 .+20;
      label = "Polyimide (Dharmadasa, 2020)",
      linewidth=8,
      color = :orange,
      clip = true)
    
    axislegend(ax1_1, labelsize=font-10, framevisible=false,position = :rb)
    fig1



## Fig3 Strain Rate Einfluss
fig3 = Figure(resolution = (1000, 600));
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

σ_exp=[]
ϵ_exp=[]
dϵdt=[]


##### Strain Rate 0.7638 %/min

BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/AP5ZK03/251006-0"
name="251006-"
df_Group=[]
df_Group_Total=[]

for i in [2]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0.01,2.5)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end

df_Group
push!(σ_exp,df_Group[1][!,:Stress]*1e6)
push!(ϵ_exp,df_Group[1][!,:Strain])
push!(dϵdt,0.007638/60)


    ##### Strain Rate 1.43 %/min

    BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/AP5ZK03/251006-0"
    name="251006-"
    df_Group=[]
    df_Group_Total=[]

    for i in [1]
            file_path = "$(BaseName)$(i).csv"
            df=import_TIRA(file_path)
            df[!,:Name]=fill("$(name)$(i)", nrow(df))
            df[!,:Stress]=df[!,:Kraft]*1e-6/As
            df[!,:Strain]=df[!,:dL_ORG]/100
            df_filtered=filter_range(df,"dL_ORG",0.01,2.5)
            push!(df_Group,df_filtered)
            push!(df_Group_Total,df_filtered)
    end

    df_Group
    push!(σ_exp,df_Group[1][!,:Stress]*1e6)
    push!(ϵ_exp,df_Group[1][!,:Strain])
    push!(dϵdt,0.0143/60)



##### Strain Rate 4.44 %/min

BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/AP5ZK03/251006-0"
name="251006-"
df_Group=[]
df_Group_Total=[]

for i in [6]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0.01,2.5)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end

df_Group
push!(σ_exp,df_Group[1][!,:Stress]*1e6)
push!(ϵ_exp,df_Group[1][!,:Strain])
push!(dϵdt,0.0444/60)


### Strain rate 27%
BaseName="D:/01 Projekt/03 Programm/Rini aktuell/Rini/development/zugversuch/AP5ZK04/251103-"
name="251103-"
df_Group=[]
df_Group_Total=[]

for i in [1]
        file_path = "$(BaseName)$(i).csv"
        df=import_TIRA(file_path)
        df[!,:Name]=fill("$(name)$(i)", nrow(df))
        df[!,:Stress]=df[!,:Kraft]*1e-6/As
        df[!,:Strain]=df[!,:dL_ORG]/100
        df_filtered=filter_range(df,"dL_ORG",0.01,2.5)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end


push!(σ_exp,df_Group[1][!,:Stress]*1e6)
push!(ϵ_exp,df_Group[1][!,:Strain])
push!(dϵdt,0.2703/60)


for i in 1:length(σ_exp)
    Makie.scatter!(ax1_3 ,ϵ_exp[i]*100,σ_exp[i]/1e6;
    markersize=14, 
    label = "v=$(round(dϵdt[i]*60*100, digits=3)) %/min")
end

 axislegend(ax1_3, labelsize=font-10, framevisible=false,position = :rb)

fig3