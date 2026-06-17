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
    
#### Log
    ## 
###############################

function simulate_stress(ParyModel,ϵ,dϵdt)

    ## Maxwell
    #k1,k2,n=[ParyModel[1],ParyModel[2],ParyModel[3]]
    #σ=dϵdt*n*(1 .- exp.(-k2.*ϵ./(n*dϵdt))) .+k1 .* ϵ 

    ## Prony
    E_inf=ParyModel[1]
    E1=ParyModel[2]
    E2=ParyModel[3]
    E3=ParyModel[4]
    τ1=ParyModel[5]
    τ2=ParyModel[6]
    τ3=ParyModel[7]
    t= ϵ./dϵdt
    E=E_inf.+E1.*exp.(−t/τ1) .+E2.*exp.(−t/τ2) .+E3.*exp.(−t/τ3) 
    
    σ=E.*ϵ
    return σ
end


# Define the predictive function

# Objective function to minimize (difference between model and experimental data)
function objective_function(ParyModel,σ_exp,ϵ_exp,dϵdt)
    total_error=0
        for i in 1:length(σ_exp)
            σ_sim    = simulate_stress(ParyModel, ϵ_exp[i],dϵdt[i])
            
            # Calculate sum of squared errors
            
            error = sum((σ_sim .- σ_exp[i]).^2)
            total_error = total_error+error
        end

    return total_error
end

# Main optimization function

function optimize_Pary(ParyModel0,σ_exp,ϵ_exp,dϵdt)
    
    max_iterations=500

    ##Maxwell
    #initial_params = [ParyModel0[1],ParyModel0[2],ParyModel0[3]]

    #lower_bounds = [1e6, 1e9,1e9]
    #upper_bounds = [ 1e15, 1e15, 1e15]

    ## Prony

    #σ=E_inf+E1*exp(−t/τ1) +E2*exp(−t/τ2) 
    lower_bounds = [1e6, 1e6,1e6,1e6,0.00001,0.00001,0.00001]
    upper_bounds = [ 1e15, 1e15, 1e15, 1e15,1e2, 1e2, 1e2]
    initial_params = [ParyModel0[1],ParyModel0[2],ParyModel0[3],ParyModel0[4],ParyModel0[5],ParyModel0[6],ParyModel0[7]]

    # Define the objective function with fixed experimental data
    obj(ParyModel_fit) = objective_function(ParyModel_fit,σ_exp,ϵ_exp,dϵdt)
    
    # Perform optimization using Nelder-Mead or L-BFGS (you can switch methods)
    result = optimize(obj, 
                     lower_bounds, 
                     upper_bounds, 
                     initial_params, 
                     Fminbox(LBFGS()),
                     Optim.Options(iterations=max_iterations, show_trace=true))
    
    # Extract optimized parameters
    optimized_params = Optim.minimizer(result)
    minimum_error = Optim.minimum(result)
    
    return optimized_params, minimum_error
end


########### Optimization algorithm ###############

function main(ParyModel0,σ_exp,ϵ_exp,dϵdt)
    
    # Run optimization
    optimized_params, error = optimize_Pary(ParyModel0,σ_exp,ϵ_exp,dϵdt)

    
    return optimized_params,error
end


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



#ParyModel0=[1e9, 1e10,10e9]
ParyModel0=[1e9, 1e9,1e9,1e9,0.1,0.1,0.1] # Prony
ParyModel_fit,error=main(ParyModel0,σ_exp,ϵ_exp,dϵdt) # Find the better adjusment for the CZM parameter
σ_sim=[]

for i in 1:length(ϵ_exp)
    σ_cum=[]
    for ii in ϵ_exp[i]
        push!(σ_cum,simulate_stress(ParyModel_fit,ii,dϵdt[i]))
    end 
    push!(σ_sim,σ_cum)
end

σ_sim

# Plot the first curve with lines and markers
fig1 = Figure(resolution = (1000, 600));
fig2 = Figure(resolution = (1000, 600));
fig3 = Figure(resolution = (1000, 600));
font=36
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




## Fig1

    Makie.scatter!(ax1_1,ϵ_exp[1]*100,σ_exp[1]/1e6;
     label = " Exp $(round(dϵdt[1]*60*100, digits=3)) [%/min]",
     marker = :utriangle,
      markersize = 14,
      color = :orange)

    Makie.lines!(ax1_1,ϵ_exp[1]*100,σ_sim[1]/1e6;
     linewidth = 4, 
     label = "Model $(round(dϵdt[1]*60*100, digits=3)) [%/min]",
     color = :black)
    
    Makie.scatter!(ax1_1,ϵ_exp[2]*100,σ_exp[2]/1e6; 
    markersize=14, 
    label = " Exp $(round(dϵdt[2]*60*100, digits=3)) [%/min]",
      marker = :rect,
      color = :red)
    Makie.lines!(ax1_1,ϵ_exp[2]*100,σ_sim[2]/1e6;
      label = "Model $(round(dϵdt[2]*60*100, digits=3)) [%/min]",
      linewidth = 4, 
      color = :black,
      linestyle = :dash)
    axislegend(ax1_1, labelsize=font-10, framevisible=false,position = :rb)


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
     label = "T22°C v=$(round(dϵdt[1]*60*100, digits=3)) [%/min]",
     color = :black)
    
    Makie.scatter!(ax1_2,ϵ_exp[2]*100,σ_exp[2]/1e6;
      label = "T50°C v=$(round(dϵdt[2]*60*100, digits=3)) [%/min]",
      markersize=14,
      color = :orange)
    
    Makie.scatter!(ax1_2,ϵ_exp[3]*100,σ_exp[3]/1e6;
    label = "T100°C v=$(round(dϵdt[3]*60*100, digits=3)) [%/min]",
    markersize=14, 
    color = :red)
    axislegend(ax1_2, labelsize=font-10, framevisible=false,position = :rb)

    fig2


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