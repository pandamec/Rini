using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Plots")
Pkg.add("Polynomials")
using CSV
using DataFrames
using Polynomials
using CairoMakie
w=5*1e-3
t=15*1e-6
l=15*1e-3

As=w*t

using Plots

function import_TIRA(file_path)
    
    df=CSV.read(file_path,DataFrame) 
    df.Zeit = parse.(Float64, replace.(df.Zeit, "," => "."))
    df.Länge = parse.(Float64, replace.(df.Länge, "," => "."))
    df.Weg = parse.(Float64, replace.(df.Weg, "," => "."))
    df.Dehnung = parse.(Float64, replace.(df.Dehnung, "," => "."))
    df.Kraft = parse.(Float64, replace.(df.Kraft, "," => "."))
    df.dL_ORG = parse.(Float64, replace.(df.dL_ORG, "," => "."))

    return df

end
# Replace with the ath to your Excel file
file_path = "test/250522-1.csv"

df=import_TIRA(file_path)
T=df[!,:Zeit]
#X=df[!,:Dehnung]
#X4=df[!,:Länge]
X=df[!,:Weg]/18
X2= df[6989:end,:Länge]/18 # According to Uwe
#X2= df[6989:end,:Dehnung] # According to Uwe

X3=df[!,:Weg]
Y=(df[6989:end,:Kraft].-df[6989,:Kraft]).*1e-6/As
Y2=df[!,:Kraft]

# Fit a 1st-degree polynomial (linear)
p_strain=fit(T[1:1000],X2[1:1000]*100*60,1)
p = fit(X2[1:1000], Y[1:1000], 1)
Xfit=X2[1:1000]
Yfit = p.(Xfit)

# Plot the first curve with lines and markers
fig = Figure(resolution = (1000,600))

font=36
ax = Axis(fig[1,1],
    xlabel = L"\epsilon (%)",
    ylabel = L"\sigma (MPa)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    #xgridstyle = :dash,        # dashed grid
    #ygridstyle = :dash
    xgridvisible = false,
    ygridvisible = false
)

Makie.scatter!(ax,X2[1:12000]*100, Y[1:12000]; markersize=12, label="Width 5 mm, strain rate 5.5 %/min",color=:orange)
#Makie.lines!(ax,Xfit*100, Yfit; linewidth=5, label= "linear fit",color=:black)
axislegend(ax, position=:rb, labelsize=font-10, framevisible=false)
fig

save("tensile.svg",fig)


# Plot the first curve with lines and markers
fig2 = Figure(resolution = (1000,600))

font=36
ax2 = Axis(fig2[1,1],
    xlabel = L"Verschiebung (mm)",
    ylabel = L"F (N)",
    xlabelsize = font,
    ylabelsize = font,
    xticklabelsize = font - 2,
    yticklabelsize = font - 2,
    #xgridstyle = :dash,        # dashed grid
    #ygridstyle = :dash
    xgridvisible = false,
    ygridvisible = false

)
p_strain=fit(T[7000:8000],X3[7000:8000]*60,1)
Makie.scatter!(ax2,X3, Y2; markersize=12, label="Width 5 mm, speed rate 1 mm/min",color=:orange)
#Makie.lines!(ax,Xfit*100, Yfit; linewidth=5, label= "linear fit",color=:black)
axislegend(ax2, position=:rb, labelsize=font-10, framevisible=false)
fig2

########




Plots.scatter(X2[1:1000]*100, Y[1:1000], label="Width 5 mm, speed mm/min", linewidth=1)
Plots.plot!(Xfit*100, Yfit, linewidth=3)


# Customize title and axes
#title!("Tensile test")
Plots.xlabel!(L"\epsilon (\%)")
Plots.ylabel!(L"\sigma (MPa)")


# Replace with the path to your Excel file
file_path = "test/250521_Zugversuch1.csv"

w=15*1e-3
t=15*1e-6
l=20*1e-3
As=w*t

df=import_TIRA(file_path)
T=df[!,:Zeit]
X=df[!,:Länge]./21
Y=df[!,:Kraft].*1e-6/As


p = fit(X[500:2500], Y[500:2500], 1)
Xfit=X[500:2500]
Yfit = p.(Xfit)


plot(X, Y, label="Width 15mm, speed 0,05mm/min", linewidth=1)
Plots.plot!(Xfit, Yfit, linewidth=1)
