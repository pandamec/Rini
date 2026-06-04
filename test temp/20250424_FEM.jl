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
file_path = "Programs/250605_Simulation_Static.xlsx"


xf=XLSX.openxlsx(file_path) 
sheet = xf[1]  # or xf["SheetName"]

# Read entire sheet content into a matrix
data = sheet[:]  # or sheet["A1:Z100"] for a specific range
    
# Convert to DataFrame; first row is assumed to be the header
header = data[1, :]
values = data[2:end, :]
df = DataFrame(values, Symbol.(header))
    
df

X=df[!,:Cracklength][1:10]
Y1=df[!,:SyParyleneMax][1:10]
Y2=df[!,:SyParyleneMax][11:20]
Y3=df[!,:SyParyleneMax][21:30]
Y4=df[!,:SyParyleneMax][31:40]
Y5=df[!,:SyParyleneMax][41:50]
# Plot the first curve with lines and markers
plot(X, Y1, label="Parylene length 200um", marker=:diamond, linewidth=2)

# Add additional curves with lines and markers
plot!(X, Y2, label="Parylene length 225um", marker=:square, linewidth=2)
plot!(X, Y3, label="Parylene length 250um", marker=:circle, linewidth=2)
plot!(X, Y4, label="Parylene length 275um", marker=:xcross, linewidth=2)
plot!(X, Y5, label="Parylene length 300um", marker=:utriangle, linewidth=2)

#plot!(x, y3, label="Curve 3", marker=:diamond, linewidth=2)

# Customize title and axes
title!("Nominal stress (static)")
xlabel!("Pre-crack length (mm)")
ylabel!("Normal stress (MPa)")


## Curvature (through displacement) influence ####

X=df[!,:Cracklength][1:10]
Y1=df[!,:SyParyleneMax][51:60]
Y2=df[!,:SyParyleneMax][61:70]
Y3=df[!,:SyParyleneMax][1:10]
Y4=df[!,:SyParyleneMax][71:80]

# Plot the first curve with lines and markers
plot(X, Y1, label="w 0.5mm", marker=:diamond, linewidth=2)

# Add additional curves with lines and markers
plot!(X, Y2, label="w 1mm", marker=:square, linewidth=2)
plot!(X, Y3, label="w 1.5mm", marker=:circle, linewidth=2)
plot!(X, Y4, label="w 2mm", marker=:xcross, linewidth=2)


# Customize title and axes
title!("Nominal Stress (static)")
xlabel!("Pre-crack length (mm)")
ylabel!("Normal stress (MPa)")



### Delamination

# Replace with the path to your Excel file
file_path = "Programs/250607_Simulation_Delamination.xlsx"


xf=XLSX.openxlsx(file_path) 
sheet = xf[1]  # or xf["SheetName"]

# Read entire sheet content into a matrix
data = sheet[:]  # or sheet["A1:Z100"] for a specific range
    
# Convert to DataFrame; first row is assumed to be the header
header = data[1, :]
values = data[2:end, :]
df2 = DataFrame(values, Symbol.(header))
    
df2

X=df2[!,:Cracklength][1:10]
Y1=df2[!,:SyParyleneMax][1:10]
Y2=df2[!,:SyParyleneMax][11:20]
Y3=df2[!,:SyParyleneMax][21:30]
Y4=df2[!,:SyParyleneMax][31:40]
Y5=df2[!,:SyParyleneMax][41:50]
# Plot the first curve with lines and markers
plot(X, Y1, label="Parylene length 200um", marker=:diamond, linewidth=2)

# Add additional curves with lines and markers
plot!(X, Y2, label="Parylene length 225um", marker=:square, linewidth=2)
plot!(X, Y3, label="Parylene length 250um", marker=:circle, linewidth=2)
plot!(X, Y4, label="Parylene length 275um", marker=:xcross, linewidth=2)
plot!(X, Y5, label="Parylene length 300um", marker=:utriangle, linewidth=2)



#plot!(x, y3, label="Curve 3", marker=:diamond, linewidth=2)

# Customize title and axes
title!("Stress under delamination")
xlabel!("Pre-crack length (mm)")
ylabel!("Normal stress (MPa)")

######### Curvature (through displacement) influence ####

X=df2[!,:Cracklength][1:10]
Y1=df2[!,:SyParyleneMax][51:60]
Y2=df2[!,:SyParyleneMax][61:70]
Y3=df2[!,:SyParyleneMax][1:10]
Y4=df2[!,:SyParyleneMax][71:80]

# Plot the first curve with lines and markers
plot(X, Y1, label="w 0.5mm", marker=:diamond, linewidth=2)

# Add additional curves with lines and markers
plot!(X, Y2, label="w 1mm", marker=:square, linewidth=2)
plot!(X, Y3, label="w 1.5mm", marker=:circle, linewidth=2)
plot!(X, Y4, label="w 2mm", marker=:xcross, linewidth=2)


# Customize title and axes
title!("Stress under delamination")
xlabel!("Pre-crack length (mm)")
ylabel!("Normal stress (MPa)")
