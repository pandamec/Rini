# 24.04.2025 Rev 0.1.0 First version
# 05.06.2025 Rev 0.1.1 Parylene lengths added 

using Pkg
Pkg.add("XLSX")
Pkg.add("DataFrames")
Pkg.add("Plots")

using XLSX
using DataFrames
using Plots
# Replace with the path to your Excel file
file_path = "test/251001_AP5S02.xlsx"

xf=XLSX.openxlsx(file_path) 
sheet = xf[1]  # or xf["SheetName"]

# Read entire sheet content into a matrix
data = sheet[:]  # or sheet["A1:Z100"] for a specific range
    
# Convert to DataFrame; first row is assumed to be the header
header = data[1, :]
values = data[2:end, :]
df = DataFrame(values, Symbol.(header))
    
df
df.ao = Float64.(df.ao)

G=(0.1*1e6*0.1*1e-3)/2
G=(0.5*1e6*0.1*1e-3)/2
G=(1*1e6*0.1*1e-3)/2
G=(5*1e6*0.1*1e-3)/2
G=(10*1e6*0.1*1e-3)/2

## CZM Parameters

δ_max=0.001e-3 #m
σ_max=[0.01 0.05 0.1 0.5 1] #MPa
G_max=σ_max*1e6*δ_max/2

df_sorted=[]
for i in σ_max

    df_sigma= df[coalesce.(df."SigmaMax" .== i, false), :]
    df_sigma[!,:σ]=df_sigma[!,:Fy]/0.1
    df_sigma[!,:δ]=df_sigma[!,:uy1]-df_sigma[!,:uy2]
    push!(df_sorted,df_sigma)


end


font = 32
using CairoMakie

fig = Figure(resolution = (1000, 600))

ax = Axis(fig[1, 1],
    xlabel = L"ao (μm)",
    ylabel = L"δ (um)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 5,
    yticklabelsize = font - 5,
    xgridvisible = false,   # kein vertikales Grid
    ygridvisible = false
    #yticks=5:5:40,
    #xticks=50:150:1000
)

#ax.xticks = 0:5000:20000


CairoMakie.scatter!(ax, df_sorted[1][!,:ao], df_sorted[1][!,:δ]*1000,markersize=15,marker=:diamond)
CairoMakie.lines!(ax, df_sorted[1][!,:ao], df_sorted[1][!,:δ]*1000, linewidth=3,label="G_max=0.005")

CairoMakie.scatter!(ax, df_sorted[2][!,:ao], df_sorted[2][!,:δ]*1000,markersize=15,marker=:diamond)
CairoMakie.lines!(ax, df_sorted[2][!,:ao], df_sorted[2][!,:δ]*1000, linewidth=3,label="G_max=0.025")

CairoMakie.scatter!(ax, df_sorted[3][!,:ao], df_sorted[3][!,:δ]*1000,markersize=15,marker=:diamond)
CairoMakie.lines!(ax, df_sorted[3][!,:ao], df_sorted[3][!,:δ]*1000, linewidth=3,label="G_max=0.05")


CairoMakie.scatter!(ax, df_sorted[4][!,:ao], df_sorted[4][!,:δ]*1000,markersize=15,marker=:diamond)
CairoMakie.lines!(ax, df_sorted[4][!,:ao], df_sorted[4][!,:δ]*1000, linewidth=3,label="G_max=0.25")


CairoMakie.scatter!(ax, df_sorted[5][!,:ao], df_sorted[5][!,:δ]*1000,markersize=15,marker=:diamond)
CairoMakie.lines!(ax, df_sorted[5][!,:ao], df_sorted[5][!,:δ]*1000, linewidth=3,label="G_max=0.5")
axislegend(ax, position = :rb) 
fig

font = 32
using CairoMakie

fig2 = Figure(resolution = (1000, 600))

ax2 = Axis(fig2[1, 1],
    xlabel = L"ao (μm)",
    ylabel = L"σ (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 5,
    yticklabelsize = font - 5,
    xgridvisible = false,   # kein vertikales Grid
    ygridvisible = false
    #yticks=5:5:40,
    #xticks=50:150:1000
)


CairoMakie.scatter!(ax2, df_sorted[1][!,:ao], df_sorted[1][!,:σ],markersize=15,marker=:diamond)
CairoMakie.lines!(ax2, df_sorted[1][!,:ao], df_sorted[1][!,:σ], linewidth=3,label="G_max=0.005")

CairoMakie.scatter!(ax2, df_sorted[2][!,:ao], df_sorted[2][!,:σ]*1,markersize=15,marker=:diamond)
CairoMakie.lines!(ax2, df_sorted[2][!,:ao], df_sorted[2][!,:σ]*1, linewidth=3,label="G_max=0.025")

CairoMakie.scatter!(ax2, df_sorted[3][!,:ao], df_sorted[3][!,:σ]*1,markersize=15,marker=:diamond)
CairoMakie.lines!(ax2, df_sorted[3][!,:ao], df_sorted[3][!,:σ]*1, linewidth=3,label="G_max=0.05")


CairoMakie.scatter!(ax2, df_sorted[4][!,:ao], df_sorted[4][!,:σ]*1,markersize=15,marker=:diamond)
CairoMakie.lines!(ax2, df_sorted[4][!,:ao], df_sorted[4][!,:σ]*1, linewidth=3,label="G_max=0.25")


CairoMakie.scatter!(ax2, df_sorted[5][!,:ao], df_sorted[5][!,:σ]*1,markersize=15,marker=:diamond)
CairoMakie.lines!(ax2, df_sorted[5][!,:ao], df_sorted[5][!,:σ]*1, linewidth=3,label="G_max=0.5")
axislegend(ax2, position = :rb) 
fig2

fig