
using CairoMakie
using DataFrames
using CSV
using Statistics

function importTest(path)

    #file = XLSX.readxlsx(path)
    df=CSV.read(path,DataFrame) 
    
    for col in names(df)
        if eltype(df[!, col]) <: AbstractString
           
                df[!, col] = parse.(Float64, replace.(df[!, col], "," => "."))
                # skip columns that cannot be converted
            
           
        end
    end

    return df
end

file_path = "development/fatigue delamination/AP5F10_Fatigue.csv"

df=importTest(file_path)
df_study=importTest(file_path)

df

delete!(df, 2)
delete!(df, 8)
delete!(df, 3)

df[!,:a]=2000 .- df[!,:lf]

df[!,:a0]=2000 .- df[!,:l0]
df_study[!,:a0]=2000 .- df_study[!,:l0]
df[!,:delta_a]= df[!,:a] -df[!,:a0]
#df[!,:N0]= zeros(length(df[!,:a0]))


df_A = df[coalesce.(df."Sample" .== "A", false), :]
df_B = df[coalesce.(df."Sample" .== "B", false), :]
df_N10 = df[coalesce.(df."N" .== 10000, false), :]
df_N20 = df[coalesce.(df."N" .== 20000, false), :]
df_N30 = df[coalesce.(df."N" .== 30000, false), :]

font = 20
fig = Figure(resolution = (1200, 800))

ax = Axis(fig[1, 1],
    xlabel = L"N",
    ylabel = L"a (um)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    xgridvisible = false,   # kein vertikales Grid
    ygridvisible = false 

)
ax.xticks = 0:5000:20000

CairoMakie.scatter!(ax, df."N", df."delta_a")

fig

save("test.pdf", fig)


# Example data
categories = ["10000", "15000", "20000"]

avg=[mean(df_N10."delta_a") , mean(df_N20."delta_a") , mean(df_N30."delta_a")]

delete!(df_study, 2)
delete!(df_study, 8)
delete!(df_study, 3)
delete!(df_study,5)
avg_crack=[mean(df_study."a0"),mean(df_N10."a") , mean(df_N20."a") , mean(df_N30."a")]

font=24
# Create figure

fig = Figure(resolution = (600, 400))

ax = Axis(fig[1, 1], 
    xlabel = L"N", ylabel = L"a(um)",
        xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font -5,
    yticklabelsize = font -5,
    xgridvisible = false,   # kein vertikales Grid
    ygridvisible = false  
)


# Bar plot
barplot!(ax, 1:length(categories), avg; color=:gray)
# Set x-tick labels to category names
errors_max=[(maximum(df_N10."delta_a") - avg[1]), (maximum(df_N15."delta_a") - avg[2]), (maximum(df_N20."delta_a") - avg[3])]
errors_min=[(minimum(df_N10."delta_a") - avg[1]), (minimum(df_N15."delta_a") - avg[2]), (minimum(df_N20."delta_a") - avg[3])]

for (i, (v, emax,emin)) in enumerate(zip(avg, errors_max,errors_min))
    lines!(ax, [i, i], [v + emin, v + emax], color=:black)
    lines!(ax, [i-0.1, i+0.1], [v + emax, v + emax], color=:black)  # top cap
    lines!(ax, [i-0.1, i+0.1], [v + emin, v + emin], color=:black)  # bottom cap
end

ax.xticks = (1:length(categories), categories)

fig

# Save figure

save("barplot.pdf", fig)


error_a0=[minimum(df."a0"), maximum(df."a0") ]
error_N10=[minimum(df_N10."a"), maximum(df_N10."a")]
error_N20=[minimum(df_N20."a"), maximum(df_N20."a")]
error_N30=[minimum(df_N30."a"), maximum(df_N30."a")]

error_max=[maximum(df_study."a0"), maximum(df_N10."a"), maximum(df_N20."a"), maximum(df_N30."a")]

error_min=[minimum(df."a0"), minimum(df_N10."a"), minimum(df_N20."a"), minimum(df_N30."a")]


font=36
fig2= Figure(resolution = (1000, 600))

ax2 = Axis(fig2[1, 1], 
    xlabel = L"N", ylabel = L"a(um)",
        xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font -10,
    yticklabelsize = font -10,
    xgridvisible = false,   # kein vertikales Grid
    ygridvisible = false  
)

categories = ["0", "10000", "20000", "30000"]

barplot!(ax2, 1:length(categories), avg_crack, color=:peachpuff,width = 0.5)

for (i, (v, emax,emin)) in enumerate(zip(avg_crack, error_max,error_min))
    lines!(ax2, [i, i], [emin,  emax], color=:black)
    lines!(ax2, [i-0.1, i+0.1], [ emax, emax], color=:black)  # top cap
    lines!(ax2, [i-0.1, i+0.1], [ emin,  emin], color=:black)  # bottom cap
end

ax2.xticks = (1:length(categories), categories)
fig2
save("my_plot.svg",fig2)
