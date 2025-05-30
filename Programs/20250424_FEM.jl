using Pkg
Pkg.add("XLSX")
Pkg.add("DataFrames")
Pkg.add("Plots")

using XLSX
using DataFrames
using Plots
# Replace with the path to your Excel file
file_path = "Programs/250424_Simulation_Static.xlsx"


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
# Plot the first curve with lines and markers
plot(X, Y1, label="Parylene length 2mm", marker=:diamond, linewidth=2)

# Add additional curves with lines and markers
plot!(X, Y2, label="Parylene length 3mm", marker=:square, linewidth=2)
#plot!(x, y3, label="Curve 3", marker=:diamond, linewidth=2)

# Customize title and axes
title!("Nominal stress (Static)")
xlabel!("Pre-crack length (mm)")
ylabel!("Normal stress (MPa)")


### Delamination

# Replace with the path to your Excel file
file_path = "Programs/250424_Simulation_Delamination.xlsx"


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
# Plot the first curve with lines and markers
plot(X, Y1, label="Parylene length 2mm", marker=:diamond, linewidth=2)

# Add additional curves with lines and markers
plot!(X, Y2, label="Parylene length 3mm", marker=:square, linewidth=2)
#plot!(x, y3, label="Curve 3", marker=:diamond, linewidth=2)

# Customize title and axes
title!("Stress under delamination")
xlabel!("Pre-crack length (mm)")
ylabel!("Normal stress (MPa)")
