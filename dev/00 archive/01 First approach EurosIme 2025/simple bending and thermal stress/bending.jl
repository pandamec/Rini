using LinearAlgebra
using QuadGK  # For numerical integration
using Plots

# Beam parameters
L = 1.0          # Length of the beam (m)
E = 2.0e11       # Young's modulus (Pa), e.g., steel
I1 = 1.0e-4      # Moment of inertia, layer 1 (m^4)
I2 = 1.2e-4      # Moment of inertia, layer 2 (m^4)
I3 = 0.8e-4      # Moment of inertia, layer 3 (m^4)
EI1 = E * I1     # Flexural rigidity, layer 1
EI2 = E * I2     # Flexural rigidity, layer 2
EI3 = E * I3     # Flexural rigidity, layer 3
h = 0.01         # Thickness between layers (m), for coupling

# Basis functions: sin(jπx/L) for simply supported beam
function phi(j, x, L)
    return sin(j * π * x / L)
end

function phi_double_prime(j, x, L)
    return -(j * π / L)^2 * sin(j * π * x / L)
end

function phi_fourth_prime(j, x, L)
    return (j * π / L)^4 * sin(j * π * x / L)
end

# Number of terms in the approximation
n_terms = 3  # Adjustable

# Assemble stiffness matrix for a single layer
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

# Coupling matrix (represents interaction between layers)
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

# Target shape (e.g., parabolic for middle layer)
function target_shape(x, L)
    return 0.01 * x * (L - x)  # Example: w_2 tries to follow this
end

# Assemble global system for three layers
function assemble_three_layer_system(n_terms, L, EI1, EI2, EI3, h)
    K1 = stiffness_matrix(n_terms, L, EI1)
    K2 = stiffness_matrix(n_terms, L, EI2)
    K3 = stiffness_matrix(n_terms, L, EI3)
    C = coupling_matrix(n_terms, L)
    
    # Coupling stiffness (simplified shear stiffness between layers, k = EI/h^3 approximation)
    k = EI2 / (h^3)  # Interlayer shear stiffness (adjustable)
    
    # Global stiffness matrix (3n x 3n)
    n = n_terms
    K_global = zeros(3n, 3n)
    K_global[1:n, 1:n] = K1 + k * C           # Layer 1: strain + coupling to layer 2
    K_global[1:n, n+1:2n] = -k * C           # Coupling to layer 2
    K_global[n+1:2n, 1:n] = -k * C           # Layer 2: coupling to layer 1
    K_global[n+1:2n, n+1:2n] = K2 + 2k * C   # Layer 2: strain + coupling to 1 and 3
    K_global[n+1:2n, 2n+1:3n] = -k * C       # Coupling to layer 3
    K_global[2n+1:3n, n+1:2n] = -k * C       # Layer 3: coupling to layer 2
    K_global[2n+1:3n, 2n+1:3n] = K3 + k * C  # Layer 3: strain + coupling to 2
    
    return K_global
end

# Load vector: only middle layer follows target shape
function load_vector(n_terms, L, target_func)
    F = zeros(3 * n_terms)
    n = n_terms
    for i = 1:n
        integrand(x) = target_func(x, L) * phi(i, x, L)
        F[n + i], _ = quadgk(integrand, 0, L, rtol=1e-6)
    end
    return F
end

# Compute deflection for a layer
function deflection(x, c, n_terms, layer_idx, L)
    n = n_terms
    start_idx = (layer_idx - 1) * n + 1
    w = 0.0
    for j = 1:n_terms
        w += c[start_idx + j - 1] * phi(j, x, L)
    end
    return w
end

# Compute distributed load q(x) = EI * w''''(x) for a layer
function distributed_load(x, c, n_terms, layer_idx, EI, L)
    n = n_terms
    start_idx = (layer_idx - 1) * n + 1
    q = 0.0
    for j = 1:n_terms
        q += c[start_idx + j - 1] * EI * phi_fourth_prime(j, x, L)
    end
    return q
end

# Main computation
K_global = assemble_three_layer_system(n_terms, L, EI1, EI2, EI3, h)
F_global = load_vector(n_terms, L, target_shape)
c_global = K_global \ F_global  # Solve for all coefficients

# Extract coefficients for each layer
n = n_terms
c1 = c_global[1:n]          # Layer 1
c2 = c_global[n+1:2n]       # Layer 2
c3 = c_global[2n+1:3n]      # Layer 3

println("Coefficients Layer 1: ", c1)
println("Coefficients Layer 2: ", c2)
println("Coefficients Layer 3: ", c3)

# Evaluate deflections and loads
x_vals = range(0, L, length=100)
w1_vals = [deflection(x, c_global, n_terms, 1, L) for x in x_vals]
w2_vals = [deflection(x, c_global, n_terms, 2, L) for x in x_vals]
w3_vals = [deflection(x, c_global, n_terms, 3, L) for x in x_vals]

# Compute interlayer loads
q12_vals = [distributed_load(x, c_global, n_terms, 1, EI1, L) - 
            distributed_load(x, c_global, n_terms, 2, EI2, L) for x in x_vals]
q23_vals = [distributed_load(x, c_global, n_terms, 2, EI2, L) - 
            distributed_load(x, c_global, n_terms, 3, EI3, L) for x in x_vals]

# Plotting
p1 = plot(x_vals, w1_vals, label="Layer 1 w1(x)", xlabel="x (m)", ylabel="Deflection (m)", title="Three-Layered Beam Deflection")
plot!(x_vals, w2_vals, label="Layer 2 w2(x)")
plot!(x_vals, w3_vals, label="Layer 3 w3(x)")

p2 = plot(x_vals, q12_vals, label="q12(x) (Layer 1-2)", xlabel="x (m)", ylabel="Load (N/m)", title="Interlayer Loads")
plot!(x_vals, q23_vals, label="q23(x) (Layer 2-3)")

display(p1)
display(p2)