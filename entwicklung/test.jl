using SymPy

x=symbols("x")  # Define symbolic variable

# Define the function to integrate
f = x^2 + 3x + 2

# Perform symbolic integration
integral_f = integrate(f, x)

println("Integral of f: ", integral_f)
