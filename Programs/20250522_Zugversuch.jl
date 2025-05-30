using Pkg
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Plots")
Pkg.add("Polynomials")
using CSV
using DataFrames
using Polynomials

w=5*1e-3
t=15*1e-6
l=15*1e-3

As=w*t


function import_TIRA(file_path)
    
    df=CSV.read(file_path,DataFrame) 
    df.Zeit = parse.(Float64, replace.(df.Zeit, "," => "."))
    df.Dehnung = parse.(Float64, replace.(df.Dehnung, "," => "."))
    df.Kraft = parse.(Float64, replace.(df.Kraft, "," => "."))

    return df

end
# Replace with the path to your Excel file
file_path = "Programs/250522-1.csv"

df=import_TIRA(file_path)
T=df[!,:Zeit]
X=df[!,:Dehnung]
Y=df[!,:Kraft].*1e-6/As
# Fit a 1st-degree polynomial (linear)

p = fit(X[7000:9000], Y[7000:9000], 1)
Xfit=X[7000:9000]
Yfit = p.(Xfit)

# Plot the first curve with lines and markers

plot(X, Y, label="Width 5mm, speed 1mm/min", linewidth=1)
plot!(Xfit, Yfit, linewidth=1)


# Customize title and axes
title!("Tensile test")
xlabel!("Strain (%)")
ylabel!("Normal Stress (MPa)")


# Replace with the path to your Excel file
file_path = "Programs/250521_Zugversuch1.csv"

w=15*1e-3
t=15*1e-6
l=20*1e-3
As=w*t

df=import_TIRA(file_path)
T=df[!,:Zeit]
X=df[!,:Dehnung]
Y=df[!,:Kraft].*1e-6/As


p = fit(X[500:2500], Y[500:2500], 1)
Xfit=X[500:2500]
Yfit = p.(Xfit)


plot(X[500:2500], Y[500:2500], label="Width 15mm, speed 0,05mm/min", linewidth=1)
plot!(Xfit, Yfit, linewidth=1)
