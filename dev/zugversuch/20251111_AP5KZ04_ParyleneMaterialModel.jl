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

function simulate_strain(ParyModel,ϵ,dϵdt)
    k1,k2,n=[ParyModel[1],ParyModel[2],ParyModel[3]]
    σ=dϵdt*n*(1 .- exp.(-k2.*ϵ./(n*dϵdt))) .+k1 .* ϵ

    return σ
end


# Define the predictive function

# Objective function to minimize (difference between model and experimental data)
function objective_function(ParyModel,σ_exp,ϵ_exp,dϵdt)

    σ_sim    = simulate_strain(ParyModel, ϵ_exp,dϵdt)
    
    # Calculate sum of squared errors
    
    error = sum((σ_sim .- σ_exp).^2)

    return error
end

# Main optimization function

function optimize_Pary(ParyModel0,σ_exp,ϵ_exp,dϵdt)
    
    max_iterations=500
    initial_params = [ParyModel0[1],ParyModel0[2],ParyModel0[3]]

    lower_bounds = [1e6, 1e9,1e9]
    upper_bounds = [ 1e15, 1e15, 1e15]

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
        df_filtered=filter_range(df,"dL_ORG",0.1,1.5)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end

df_Group
σ_exp=df_Group[1][!,:Stress]*1e6
ϵ_exp=df_Group[1][!,:Strain]
dϵdt=0.24/60

ParyModel0=[1e9, 1e10,10e9]
ParyModel_fit1,error1=main(ParyModel0,σ_exp,ϵ_exp,dϵdt) # Find the better adjusment for the CZM parameter

σ_sim=[]
for i in ϵ_exp

push!(σ_sim,simulate_strain(ParyModel_fit1,i,dϵdt))

end

σ_sim

# Plot the first curve with lines and markers
fig = Figure(resolution = (1000, 600))
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

Makie.scatter!(ax1,ϵ_exp*100,σ_exp/1e6; markersize=12, label = "Experiment")
            
Makie.scatter!(ax1,ϵ_exp*100,σ_sim/1e6; markersize=12, label = "Maxwell Modell")
axislegend(ax1, labelsize=font-10, framevisible=false)
fig



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
        df_filtered=filter_range(df,"dL_ORG",0.1,1.5)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end

df_Group
σ_exp=df_Group[1][!,:Stress]*1e6
ϵ_exp=df_Group[1][!,:Strain]
dϵdt=0.0143/60

ParyModel0=[1e9, 1e10,10e9]
ParyModel_fit2=main(ParyModel0,σ_exp,ϵ_exp,dϵdt) # Find the better adjusment for the CZM parameter

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
        df_filtered=filter_range(df,"dL_ORG",0.1,1.5)
        push!(df_Group,df_filtered)
        push!(df_Group_Total,df_filtered)
end

df_Group
σ_exp=df_Group[1][!,:Stress]*1e6
ϵ_exp=df_Group[1][!,:Strain]
dϵdt=0.0444/60

ParyModel0=[1e9, 1e10,10e9]
ParyModel_fit3=main(ParyModel0,σ_exp,ϵ_exp,dϵdt) # Find the better adjusment for the CZM parameter

####