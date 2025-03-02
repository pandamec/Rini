using CSV
using DataFrames
using Plots

function parse_float(x)
    # Replace comma with period and convert to Float64, handle scientific notation
    return parse(Float64, replace(string(x), "," => "."))
end



file_path = "D:/01 Projekt/03 Programm/Rini/Rini/datei/26022025_Sim.csv"  # Update with actual file path
# Read the CSV file into a DataFrame
# Replace "your_file.csv" with the actual path to your CSV file
df = CSV.read(file_path, DataFrame; delim=";")

numeric_columns = filter(col -> col ∉ ["Name", "Note"], names(df))
for col in numeric_columns
    df[!, col] = parse_float.(df[!, col])
end


filtered_df = filter(row -> row[Symbol("P28 - PaLange")] == 2, df)

# Step 3: Create a scatter plot
# Example: Plot P3_Cracklength (x-axis) vs P7_YAxisNormalStressKlebstoffEndTimeMaximum (y-axis)


x = filtered_df[!, Symbol("P3 - Cracklength")]
y = filtered_df[!, Symbol("P8 - Y Axis - Normal Stress - Parylene - End Time  Maximum")]

#scatter(x, y,
 #   xlabel = "Crack Length(mm)",
 #   ylabel = "Maximum stress at the interface (MPa)",
 #   label = "Lange2",
 #   markersize = 5,
 #   legend = :topright
#)

filtered2_df = filter(row -> row[Symbol("P28 - PaLange")] == 3, df)

x2 = filtered2_df[!, Symbol("P3 - Cracklength")]
y2 = filtered2_df[!, Symbol("P8 - Y Axis - Normal Stress - Parylene - End Time  Maximum")]

scatter!(x2, y2,
    xlabel = "Crack length(mm)",
    ylabel = "Maximum stress at the interface (MPa)",
    label = "Lange 3",
    markersize = 5,
    legend = :topright
)

filtered3_df = filter(row -> row[Symbol("P28 - PaLange")] == 3.5, df)

x3 = filtered3_df[!, Symbol("P3 - Cracklength")]
y3 = filtered3_df[!, Symbol("P8 - Y Axis - Normal Stress - Parylene - End Time  Maximum")]

scatter!(x3, y3,
    xlabel = "Crack length(mm)",
    ylabel = "Maximum stress at the interface (MPa)",
    label = "Lange 3.5",
    markersize = 5,
    legend = :topright
)

filtered4_df = filter(row -> row[Symbol("P28 - PaLange")] == 4, df)

x4 = filtered4_df[!, Symbol("P3 - Cracklength")]
y4 = filtered4_df[!, Symbol("P8 - Y Axis - Normal Stress - Parylene - End Time  Maximum")]

scatter!(x4, y4,
    xlabel = "Crack length(mm)",
    ylabel = "Maximum stress at the interface (MPa)",
    label = "Lange 4",
    markersize = 5,
    legend = :topright
)
