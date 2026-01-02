using PyPlot  
using Plots# Requires Pkg.add("PyPlot")

# Material properties
struct Material
    E::Float64    # Young's modulus (Pa)
    ν::Float64    # Poisson's ratio
    Gc::Float64   # Critical energy release rate (J/m²)
end

# Cohesive zone parameters
struct CohesiveLaw
    σ_max::Float64  # Maximum cohesive strength (Pa)
    δ_c::Float64    # Critical separation (m)
    Gc::Float64     # Critical energy release rate (J/m²)
end

# DCB geometry and loading
struct DCB
    L::Float64      # Total length (m)
    a0::Float64     # Initial crack length (m)
    h::Float64      # Half-thickness of each layer (m)
    b::Float64      # Width (m)
    P_max::Float64  # Maximum applied load per cycle (N)
end

# Define materials
Si = Material(170e9, 0.28, 50.0)         # Silicon
Parylene = Material(2.8e9, 0.4, 10.0)    # Parylene
interface = CohesiveLaw(50e6, 1e-6, 20.0)  # Si-Parylene interface
beam = DCB(0.1, 0.02, 0.001, 0.01, 10.0)  # DCB geometry

# Fatigue parameters (Paris law)
C = 1e-5  # Paris law coefficient (m/cycle/(J/m²)^m)
m = 2.5   # Paris law exponent

# Energy release rate for DCB (Mode I)
function energy_release_rate(P, a, h, E, b)
    I = b * h^3 / 12
    G = (3 * P^2 * a^2) / (2 * b * E * I)
    return G
end

# Cohesive zone length (Turon’s estimate)
function cohesive_zone_length(E, Gc, σ_max)
    # Turon et al. (2007): L_cz = M * E * Gc / σ_max^2, M ≈ 0.5-1 (typically 1 for simplicity)
    return E * Gc / σ_max^2
end

# Effective modulus
function effective_modulus(E1, E2)
    return (E1 + E2) / 2
end

# Fatigue simulation with Turon’s damage evolution
function simulate_fatigue_delamination(dcb, mat1, mat2, cohesive, cycles)
    a = dcb.a0
    crack_lengths = [a]
    G_values = Float64[]
    D_values = [0.0]  # Damage variable (0 = intact, 1 = fully damaged)
    
    E_eff = effective_modulus(mat1.E, mat2.E)
    L_cz = cohesive_zone_length(E_eff, cohesive.Gc, cohesive.σ_max)
    println("Cohesive zone length: $L_cz m")
    
    for n in 1:cycles
        G = energy_release_rate(dcb.P_max, a, dcb.h, E_eff, dcb.b)
        push!(G_values, G)
        
        # Fatigue damage rate (Turon’s approach)
        da_dN = C * G^m  # Paris law crack growth rate
        ΔD = (da_dN / L_cz) * (G / cohesive.Gc)  # Damage increment per cycle
        D = min(D_values[end] + ΔD, 1.0)  # Cap at 1
        push!(D_values, D)
        
        # Update crack length when fully damaged
        if D >= 1.0
            a += L_cz  # Advance crack by cohesive zone length
            push!(crack_lengths, a)
            push!(D_values, 0.0)  # Reset damage for new cohesive zone
            push!(G_values, G)   # Repeat G for continuity
        else
            push!(crack_lengths, a)  # No advance yet
        end
        
        if a >= dcb.L
            println("Crack reached beam end at cycle $n")
            break
        end
    end
    
    return crack_lengths, G_values, D_values
end

# Run simulation
max_cycles = 10000
println("Running simulation for $max_cycles cycles...")
crack_lengths, G_values, D_values = simulate_fatigue_delamination(beam, Si, Parylene, interface, max_cycles)

# Plot results
figure()

subplot(1, 3, 1)
plot(0:length(crack_lengths)-1, crack_lengths .* 1000, linewidth=2)
xlabel("Cycles")
ylabel("Crack Length (mm)")
title("Delamination Growth")

subplot(1, 3, 2)
plot(0:length(G_values)-1, G_values, linewidth=2)
xlabel("Cycles")
ylabel("Energy Release Rate (J/m²)")
title("Energy Release Rate")

subplot(1, 3, 3)
plot(0:length(D_values)-1, D_values, linewidth=2)
xlabel("Cycles")
ylabel("Damage Variable")
title("Damage Evolution")

