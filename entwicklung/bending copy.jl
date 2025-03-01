using LinearAlgebra
using QuadGK
using Plots

# Beam parameters
L = 0.002          # Length (m)
E1 = 1.2e11       # Young's modulus (Pa)
E2 = 2.8e9       # Young's modulus (Pa)
E3 = 1e9       # Young's modulus (Pa)
I1 = 0.675e-3*0.002      # Moment of inertia, layer 1 (m^4)
I2 = 0.01e-3*0.02      # Moment of inertia, layer 2 (m^4)
I3 = 0.1e-3*0.02      # Moment of inertia, layer 3 (m^4)
EI1 = E1 * I1
EI2 = E2 * I2
EI3 = E3 * I3
h = 0.675e-3         # Thickness between layers (m)
k = EI2 / (h^3)  # Interlayer shear stiffness (N/m^2)



# Beam parameters
L = 0.002          # Length (m)
b = 0.1          # Width (m)
h1 = 0.675e-3       # Thickness, layer 1 (m)
h2 = 0.01e-3        # Thickness, layer 2 (m)
h3 = 0.1e-3        # Thickness, layer 3 (m)
E1 = 1.2e11     # Young's modulus, layer 1 (Pa)
E2 = 2.8e9     # Young's modulus, layer 2 (Pa)
E3 = 1e9      # Young's modulus, layer 3 (Pa)

# Moment of inertia
I1 = b * h1^3 / 12
I2 = b * h2^3 / 12
I3 = b * h3^3 / 12

# Flexural rigidity
EI1 = E1 * I1
EI2 = E2 * I2
EI3 = E3 * I3

# Interlayer shear stiffness for each interface
E_avg12 = (E1 + E2) / 2
E_avg23 = (E2 + E3) / 2
h_avg12 = (h1 + h2) / 2
h_avg23 = (h2 + h3) / 2
k12 = (E_avg12 * h_avg12) / (h1 + h2)  # Between layer 1 and 2
k23 = (E_avg23 * h_avg23) / (h2 + h3)  # Between layer 2 and 3

# Basis functions: Simply supported
function phi(j, x, L)
    return sin(j * π * x / L)
end

function phi_double_prime(j, x, L)
    return -(j * π / L)^2 * sin(j * π * x / L)
end

function phi_fourth_prime(j, x, L)
    return (j * π / L)^4 * sin(j * π * x / L)
end

# Number of terms
n_terms = 3

# Stiffness matrix for a layer
function stiffness_matrix(n_terms, L, EI)
    K = zeros(n_terms, n_terms)
    for i = 1:n_terms
        for j = 1:n_terms
            integrand(x) = EI * phi_double_prime(i, x, L) * phi_double_prime(j, x, L)
            K[i,j], _ = quadgk(integrand, 0, L, rtol=1e-6)
        end
    end
    return K
end

# Coupling matrix
function coupling_matrix(n_terms, L)
    C = zeros(n_terms, n_terms)
    for i = 1:n_terms
        for j = 1:n_terms
            integrand(x) = phi(i, x, L) * phi(j, x, L)
            C[i,j], _ = quadgk(integrand, 0, L, rtol=1e-6)
        end
    end
    return C
end

# Assemble global system with distinct k12 and k23
function assemble_system(n_terms, L, EI1, EI2, EI3, k12, k23)
    K1 = stiffness_matrix(n_terms, L, EI1)
    K2 = stiffness_matrix(n_terms, L, EI2)
    K3 = stiffness_matrix(n_terms, L, EI3)
    C = coupling_matrix(n_terms, L)
    
    n = n_terms
    K_global = zeros(3n, 3n)
    # Layer 1: strain + coupling to layer 2
    K_global[1:n, 1:n] = K1 + k12 * C
    K_global[1:n, n+1:2n] = -k12 * C
    # Layer 2: strain + coupling to 1 and 3
    K_global[n+1:2n, 1:n] = -k12 * C
    K_global[n+1:2n, n+1:2n] = K2 + (k12 + k23) * C
    K_global[n+1:2n, 2n+1:3n] = -k23 * C
    # Layer 3: strain + coupling to 2
    K_global[2n+1:3n, n+1:2n] = -k23 * C
    K_global[2n+1:3n, 2n+1:3n] = K3 + k23 * C
    
    return K_global
end

# Target shape for layer 3
function target_shape(x, L)
    return 0.01 * x * (L - x)  # Parabolic deflection
end

# Load vector: target on layer 3
function load_vector(n_terms, L, target_func)
    F = zeros(3 * n_terms)
    n = n_terms
    for i = 1:n
        integrand(x) = target_func(x, L) * phi(i, x, L)
        F[2n + i], _ = quadgk(integrand, 0, L, rtol=1e-6)
    end
    return F
end

# Deflection for a layer
function deflection(x, c, n_terms, layer_idx, L)
    n = n_terms
    start_idx = (layer_idx - 1) * n + 1
    w = 0.0
    for j = 1:n_terms
        w += c[start_idx + j - 1] * phi(j, x, L)
    end
    return w
end

# Fourth derivative
function fourth_derivative(x, c, n_terms, layer_idx, L)
    n = n_terms
    start_idx = (layer_idx - 1) * n + 1
    w4 = 0.0
    for j = 1:n_terms
        w4 += c[start_idx + j - 1] * phi_fourth_prime(j, x, L)
    end
    return w4
end

# Main computation
K_global = assemble_system(n_terms, L, EI1, EI2, EI3, k12, k23)
F_global = load_vector(n_terms, L, target_shape)
c_global = K_global \ F_global

# Extract coefficients
n = n_terms
c1 = c_global[1:n]
c2 = c_global[n+1:2n]
c3 = c_global[2n+1:3n]

# Evaluate deflections
x_vals = range(0, L, length=100)
w1_vals = [deflection(x, c_global, n_terms, 1, L) for x in x_vals]
w2_vals = [deflection(x, c_global, n_terms, 2, L) for x in x_vals]
w3_vals = [deflection(x, c_global, n_terms, 3, L) for x in x_vals]

# Interlayer loads
q12_vals = [k12 * (w1 - w2) for (w1, w2) in zip(w1_vals, w2_vals)]
q23_vals = [k23 * (w2 - w3) for (w2, w3) in zip(w2_vals, w3_vals)]

# External load at bottom of layer 3
q3_vals = [EI3 * fourth_derivative(x, c_global, n_terms, 3, L) + 
           k23 * (deflection(x, c_global, n_terms, 2, L) - 
                  deflection(x, c_global, n_terms, 3, L)) for x in x_vals]

# Plotting
p1 = plot(x_vals, w1_vals, label="Layer 1 w1(x)", xlabel="x (m)", ylabel="Deflection (m)", title="Deflections")
plot!(x_vals, w2_vals, label="Layer 2 w2(x)")
plot!(x_vals, w3_vals, label="Layer 3 w3(x)")

p2 = plot(x_vals, q12_vals, label="q12(x) (Layer 1-2)", xlabel="x (m)", ylabel="Load (N/m)", title="Loads")
plot!(x_vals, q23_vals, label="q23(x) (Layer 2-3)")
plot!(x_vals, q3_vals, label="q3(x) (Bottom Load)")

display(p1)
display(p2)

# Print layer properties
println("Layer 1: E1 = $E1 Pa, h1 = $h1 m, I1 = $I1 m^4, EI1 = $EI1 N·m²")
println("Layer 2: E2 = $E2 Pa, h2 = $h2 m, I2 = $I2 m^4, EI2 = $EI2 N·m²")
println("Layer 3: E3 = $E3 Pa, h3 = $h3 m, I3 = $I3 m^4, EI3 = $EI3 N·m²")
println("k12 = $k12 N/m² (Layer 1-2), k23 = $k23 N/m² (Layer 2-3)")