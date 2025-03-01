using Plots
using LinearAlgebra

# Material and geometry properties
struct Material
    E::Float64      # Young's modulus (Pa)
    ν::Float64      # Poisson's ratio
    h::Float64      # Thickness of each layer (m)
    L::Float64      # Length of the beam (m)
end

struct CohesiveProperties
    σ_max::Float64  # Maximum normal strength (Pa)
    τ_max::Float64  # Maximum shear strength (Pa)
    δ_n_c::Float64  # Critical normal separation (m)
    δ_s_c::Float64  # Critical shear separation (m)
    K_n::Float64    # Normal stiffness (N/m³)
    K_s::Float64    # Shear stiffness (N/m³)
    η::Float64      # Fatigue degradation parameter (per cycle)
end

# Define properties
beam = Material(70e9, 0.3, 0.01, 0.1)  # Aluminum-like, 10 cm long, 1 cm thick
cohesive = CohesiveProperties(50e6, 30e6, 0.001, 0.002, 50e6/0.001, 30e6/0.002, 1e-5)

# Discretization
n_elements = 100
dx = beam.L / (n_elements - 1)  # Adjusted to have n_elements points
x = range(0, beam.L, length=n_elements)  # Consistent size

# Beam deflection and shear displacement (SLB configuration)
function beam_displacements(x, P, beam, a)
    I = beam.h^3 / 12  # Moment of inertia for one layer
    w = zeros(length(x))  # Vertical deflection (Mode I)
    u = zeros(length(x))  # Horizontal shear displacement (Mode II)
    
    for i in 1:length(x)
        if x[i] <= a  # Before crack tip
            w[i] = P * x[i]^2 * (3beam.L - x[i]) / (6 * beam.E * I)
            u[i] = P * x[i]^2 / (2 * beam.E * I) * beam.h / 2
        else  # After crack tip
            w[i] = P * beam.L^2 * (3x[i] - beam.L) / (6 * beam.E * I)
            u[i] = P * beam.L^2 / (2 * beam.E * I) * beam.h / 2
        end
    end
    return w, u
end

# Mixed-mode cohesive law (bilinear)
function cohesive_traction(δ_n, δ_s, cohesive, damage)
    σ_max = cohesive.σ_max * (1 - damage)
    τ_max = cohesive.τ_max * (1 - damage)
    K_n = cohesive.K_n * (1 - damage)
    K_s = cohesive.K_s * (1 - damage)
    δ_n_c = cohesive.δ_n_c
    δ_s_c = cohesive.δ_s_c
    
    # Effective separation
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(δ_n_c^2 + δ_s_c^2)
    
    if δ_eff <= δ_eff_c * σ_max / K_n
        σ = K_n * δ_n
        τ = K_s * δ_s
    elseif δ_eff < δ_eff_c
        λ = (δ_eff_c - δ_eff) / (δ_eff_c - σ_max / K_n)
        σ = σ_max * λ * (δ_n / δ_eff)
        τ = τ_max * λ * (δ_s / δ_eff)
    else
        σ = 0.0
        τ = 0.0
    end
    return σ, τ
end

# Fatigue damage evolution (mixed-mode)
function update_damage(damage, δ_n, δ_s, cohesive, cycles)
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(cohesive.δ_n_c^2 + cohesive.δ_s_c^2)
    dD_dn = cohesive.η * (δ_eff / δ_eff_c)^2
    return min(damage + dD_dn * cycles, 1.0)
end

# SLB simulation
function simulate_SLB(P, cycles, beam, cohesive)
    a = beam.L / 4  # Initial crack length
    damage = zeros(n_elements)
    δ_n_history = zeros(n_elements)
    δ_s_history = zeros(n_elements)
    σ_history = zeros(n_elements)
    τ_history = zeros(n_elements)
    
    for n in 1:cycles
        # Calculate displacements
        w, u = beam_displacements(x, P, beam, a)
        δ_n = w  # Normal separation
        δ_s = u  # Shear separation
        
        # Update tractions and damage
        for i in 1:n_elements
            σ, τ = cohesive_traction(δ_n[i], δ_s[i], cohesive, damage[i])
            σ_history[i] = σ
            τ_history[i] = τ
            if x[i] > a
                damage[i] = update_damage(damage[i], δ_n[i], δ_s[i], cohesive, 1)
            end
        end
        
        # Update crack length
        for i in 1:n_elements
            if damage[i] >= 1.0 && x[i] > a
                a = x[i]
                break
            end
        end
        
        # Store history
        δ_n_history .= δ_n
        δ_s_history .= δ_s
    end
    
    return x, δ_n_history, δ_s_history, σ_history, τ_history, damage, a
end

# Run simulation
P = 100.0  # Applied load (N)
cycles = 1000  # Number of fatigue cycles
x, δ_n, δ_s, σ, τ, damage, crack_length = simulate_SLB(P, cycles, beam, cohesive)

# Plot results
p1 = plot(x, δ_n, label="Normal Separation (m)", xlabel="Position (m)")
p2 = plot(x, δ_s, label="Shear Separation (m)", xlabel="Position (m)")
p3 = plot(x, σ, label="Normal Traction (Pa)", xlabel="Position (m)")
p4 = plot(x, τ, label="Shear Traction (Pa)", xlabel="Position (m)")
p5 = plot(x, damage, label="Damage", xlabel="Position (m)")
plot(p1, p2, p3, p4, p5, layout=(5,1), title="Cohesive Zone Model - SLB Fatigue")

##################
using Plots
using LinearAlgebra

# Material and geometry properties
struct Material
    E::Float64      # Young's modulus (Pa)
    ν::Float64      # Poisson's ratio
    h::Float64      # Thickness of each layer (m)
    L::Float64      # Length of the beam (m)
end

struct CohesiveProperties
    σ_max::Float64  # Maximum normal strength (Pa)
    τ_max::Float64  # Maximum shear strength (Pa)
    δ_n_c::Float64  # Critical normal separation (m)
    δ_s_c::Float64  # Critical shear separation (m)
    K_n::Float64    # Normal stiffness (N/m³)
    K_s::Float64    # Shear stiffness (N/m³)
    η::Float64      # Fatigue degradation parameter (per cycle)
end

# Define properties
beam = Material(70e9, 0.3, 0.01, 0.1)  # Aluminum-like, 10 cm long, 1 cm thick
cohesive = CohesiveProperties(50e6, 30e6, 0.001, 0.002, 50e6/0.001, 30e6/0.002, 1e-5)

# Discretization
n_elements = 100
x = range(0, beam.L, length=n_elements)  # Consistent size

# Beam deflection and shear displacement (SLB configuration)
function beam_displacements(x, P, beam, a)
    I = beam.h^3 / 12  # Moment of inertia for one layer
    w = zeros(length(x))  # Vertical deflection (Mode I)
    u = zeros(length(x))  # Horizontal shear displacement (Mode II)
    
    for i in 1:length(x)
        if x[i] <= a  # Before crack tip
            w[i] = P * x[i]^2 * (3beam.L - x[i]) / (6 * beam.E * I)
            u[i] = P * x[i]^2 / (2 * beam.E * I) * beam.h / 2
        else  # After crack tip
            w[i] = P * beam.L^2 * (3x[i] - beam.L) / (6 * beam.E * I)
            u[i] = P * beam.L^2 / (2 * beam.E * I) * beam.h / 2
        end
    end
    return w, u
end

# Mixed-mode cohesive law (bilinear)
function cohesive_traction(δ_n, δ_s, cohesive, damage)
    σ_max = cohesive.σ_max * (1 - damage)
    τ_max = cohesive.τ_max * (1 - damage)
    K_n = cohesive.K_n * (1 - damage)
    K_s = cohesive.K_s * (1 - damage)
    δ_n_c = cohesive.δ_n_c
    δ_s_c = cohesive.δ_s_c
    
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(δ_n_c^2 + δ_s_c^2)
    
    if δ_eff <= δ_eff_c * σ_max / K_n
        σ = K_n * δ_n
        τ = K_s * δ_s
    elseif δ_eff < δ_eff_c
        λ = (δ_eff_c - δ_eff) / (δ_eff_c - σ_max / K_n)
        σ = σ_max * λ * (δ_n / δ_eff)
        τ = τ_max * λ * (δ_s / δ_eff)
    else
        σ = 0.0
        τ = 0.0
    end
    return σ, τ
end

# Fatigue damage evolution (mixed-mode)
function update_damage(damage, δ_n, δ_s, cohesive, cycles)
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(cohesive.δ_n_c^2 + cohesive.δ_s_c^2)
    dD_dn = cohesive.η * (δ_eff / δ_eff_c)^2
    return min(damage + dD_dn * cycles, 1.0)
end

# SLB simulation with cycle tracking
function simulate_SLB(P, cycles, beam, cohesive)
    a = beam.L / 4  # Initial crack length
    damage = zeros(n_elements)
    
    # Arrays to store history at specific cycles
    cycle_points = [0, 200, 400, 600, 800, 1000]  # Track these cycles
    δ_n_history = zeros(n_elements, length(cycle_points))
    δ_s_history = zeros(n_elements, length(cycle_points))
    σ_history = zeros(n_elements, length(cycle_points))
    τ_history = zeros(n_elements, length(cycle_points))
    damage_history = zeros(n_elements, length(cycle_points))
    crack_lengths = zeros(length(cycle_points))
    
    history_idx = 1
    for n in 0:cycles
        # Calculate displacements
        w, u = beam_displacements(x, P, beam, a)
        δ_n = w
        δ_s = u
        
        # Update tractions and damage
        σ = zeros(n_elements)
        τ = zeros(n_elements)
        for i in 1:n_elements
            σ[i], τ[i] = cohesive_traction(δ_n[i], δ_s[i], cohesive, damage[i])
            if x[i] > a
                damage[i] = update_damage(damage[i], δ_n[i], δ_s[i], cohesive, 1)
            end
        end
        
        # Update crack length
        for i in 1:n_elements
            if damage[i] >= 1.0 && x[i] > a
                a = x[i]
                break
            end
        end
        
        # Store data at specified cycle points
        if n in cycle_points
            δ_n_history[:, history_idx] .= δ_n
            δ_s_history[:, history_idx] .= δ_s
            σ_history[:, history_idx] .= σ
            τ_history[:, history_idx] .= τ
            damage_history[:, history_idx] .= damage
            crack_lengths[history_idx] = a
            history_idx += 1
        end
    end
    
    return x, δ_n_history, δ_s_history, σ_history, τ_history, damage_history, cycle_points, crack_lengths
end

# Run simulation
P = 1000.0  # Increased load for visible effects (N)
cycles = 1000
x, δ_n_hist, δ_s_hist, σ_hist, τ_hist, damage_hist, cycle_points, crack_lengths = simulate_SLB(P, cycles, beam, cohesive)

# Plot results for different cycles
p1 = plot(xlabel="Position (m)", ylabel="Normal Separation (m)", title="Normal Separation")
p2 = plot(xlabel="Position (m)", ylabel="Shear Separation (m)", title="Shear Separation")
p3 = plot(xlabel="Position (m)", ylabel="Normal Traction (Pa)", title="Normal Traction")
p4 = plot(xlabel="Position (m)", ylabel="Shear Traction (Pa)", title="Shear Traction")
p5 = plot(xlabel="Position (m)", ylabel="Damage", title="Damage")

for (i, n) in enumerate(cycle_points)
    plot!(p1, x, δ_n_hist[:, i], label="Cycle $n")
    plot!(p2, x, δ_s_hist[:, i], label="Cycle $n")
    plot!(p3, x, σ_hist[:, i], label="Cycle $n")
    plot!(p4, x, τ_hist[:, i], label="Cycle $n")
    plot!(p5, x, damage_hist[:, i], label="Cycle $n")
end

plot(p1, p2, p3, p4, p5, layout=(5,1), size=(800, 1000), legend=:topleft)


#########


using Plots
using LinearAlgebra

# Material and geometry properties
struct Material
    E::Float64      # Young's modulus (Pa)
    ν::Float64      # Poisson's ratio
    h::Float64      # Thickness of each layer (m)
    L::Float64      # Length of the beam (m)
end

struct CohesiveProperties
    σ_max::Float64  # Maximum normal strength (Pa)
    τ_max::Float64  # Maximum shear strength (Pa)
    δ_n_c::Float64  # Critical normal separation (m)
    δ_s_c::Float64  # Critical shear separation (m)
    K_n::Float64    # Normal stiffness (N/m³)
    K_s::Float64    # Shear stiffness (N/m³)
    η::Float64      # Fatigue degradation parameter (per cycle)
end

# Define properties (increased η for faster fatigue)
beam = Material(70e9, 0.3, 0.01, 0.1)
cohesive = CohesiveProperties(50e6, 30e6, 0.001, 0.002, 50e6/0.001, 30e6/0.002, 1e-3)  # η = 1e-3

# Discretization
n_elements = 100
x = range(0, beam.L, length=n_elements)

# Beam deflection and shear displacement (SLB configuration)
function beam_displacements(x, P, beam, a)
    I = beam.h^3 / 12
    w = zeros(length(x))
    u = zeros(length(x))
    
    for i in 1:length(x)
        if x[i] <= a
            w[i] = P * x[i]^2 * (3beam.L - x[i]) / (6 * beam.E * I)
            u[i] = P * x[i]^2 / (2 * beam.E * I) * beam.h / 2
        else
            w[i] = P * beam.L^2 * (3x[i] - beam.L) / (6 * beam.E * I)
            u[i] = P * beam.L^2 / (2 * beam.E * I) * beam.h / 2
        end
    end
    return w, u
end

# Mixed-mode cohesive law (bilinear)
function cohesive_traction(δ_n, δ_s, cohesive, damage)
    σ_max = cohesive.σ_max * (1 - damage)
    τ_max = cohesive.τ_max * (1 - damage)
    K_n = cohesive.K_n * (1 - damage)
    K_s = cohesive.K_s * (1 - damage)
    δ_n_c = cohesive.δ_n_c
    δ_s_c = cohesive.δ_s_c
    
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(δ_n_c^2 + δ_s_c^2)
    
    if δ_eff <= δ_eff_c * σ_max / K_n
        σ = K_n * δ_n
        τ = K_s * δ_s
    elseif δ_eff < δ_eff_c
        λ = (δ_eff_c - δ_eff) / (δ_eff_c - σ_max / K_n)
        σ = σ_max * λ * (δ_n / δ_eff)
        τ = τ_max * λ * (δ_s / δ_eff)
    else
        σ = 0.0
        τ = 0.0
    end
    return σ, τ
end

# Fatigue damage evolution
function update_damage(damage, δ_n, δ_s, cohesive, cycles)
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(cohesive.δ_n_c^2 + cohesive.δ_s_c^2)
    dD_dn = cohesive.η * (δ_eff / δ_eff_c)^2
    return min(damage + dD_dn * cycles, 1.0)
end

# SLB simulation with cycle tracking
function simulate_SLB(P, cycles, beam, cohesive)
    a = beam.L / 4  # Initial crack length
    damage = zeros(n_elements)
    
    cycle_points = [0, 200, 400, 600, 800, 1000]
    δ_n_history = zeros(n_elements, length(cycle_points))
    δ_s_history = zeros(n_elements, length(cycle_points))
    σ_history = zeros(n_elements, length(cycle_points))
    τ_history = zeros(n_elements, length(cycle_points))
    damage_history = zeros(n_elements, length(cycle_points))
    crack_lengths = zeros(length(cycle_points))
    
    history_idx = 1
    for n in 0:cycles
        w, u = beam_displacements(x, P, beam, a)
        δ_n = w
        δ_s = u
        
        σ = zeros(n_elements)
        τ = zeros(n_elements)
        for i in 1:n_elements
            σ[i], τ[i] = cohesive_traction(δ_n[i], δ_s[i], cohesive, damage[i])
            if x[i] > a
                damage[i] = update_damage(damage[i], δ_n[i], δ_s[i], cohesive, 1)
            end
        end
        
        for i in 1:n_elements
            if damage[i] >= 1.0 && x[i] > a
                a = x[i]
                break
            end
        end
        
        if n in cycle_points
            δ_n_history[:, history_idx] .= δ_n
            δ_s_history[:, history_idx] .= δ_s
            σ_history[:, history_idx] .= σ
            τ_history[:, history_idx] .= τ
            damage_history[:, history_idx] .= damage
            crack_lengths[history_idx] = a
            println("Cycle $n: Crack length = $a m")  # Debug output
            history_idx += 1
        end
    end
    
    return x, δ_n_history, δ_s_history, σ_history, τ_history, damage_history, cycle_points, crack_lengths
end

# Run simulation
P = 5000.0  # Increased load for visible crack growth
cycles = 1000
x, δ_n_hist, δ_s_hist, σ_hist, τ_hist, damage_hist, cycle_points, crack_lengths = simulate_SLB(P, cycles, beam, cohesive)

# Plot results with zoom near crack tip
xlims_zoom = (0.02, 0.05)  # Focus on initial crack region
p1 = plot(xlabel="Position (m)", ylabel="Normal Separation (m)", title="Normal Separation", xlims=xlims_zoom)
p2 = plot(xlabel="Position (m)", ylabel="Shear Separation (m)", title="Shear Separation", xlims=xlims_zoom)
p3 = plot(xlabel="Position (m)", ylabel="Normal Traction (Pa)", title="Normal Traction", xlims=xlims_zoom)
p4 = plot(xlabel="Position (m)", ylabel="Shear Traction (Pa)", title="Shear Traction", xlims=xlims_zoom)
p5 = plot(xlabel="Position (m)", ylabel="Damage", title="Damage", xlims=xlims_zoom)

for (i, n) in enumerate(cycle_points)
    plot!(p1, x, δ_n_hist[:, i], label="Cycle $n")
    plot!(p2, x, δ_s_hist[:, i], label="Cycle $n")
    plot!(p3, x, σ_hist[:, i], label="Cycle $n")
    plot!(p4, x, τ_hist[:, i], label="Cycle $n")
    plot!(p5, x, damage_hist[:, i], label="Cycle $n")
end

plot(p1, p2, p3, p4, p5, layout=(5,1), size=(800, 1000), legend=:topleft)



######

using Plots
using LinearAlgebra

# Material and geometry properties
struct Material
    E::Float64      # Young's modulus (Pa)
    ν::Float64      # Poisson's ratio
    h::Float64      # Thickness of each layer (m)
    L::Float64      # Length of the beam (m)
end

struct CohesiveProperties
    σ_max::Float64  # Maximum normal strength (Pa)
    τ_max::Float64  # Maximum shear strength (Pa)
    δ_n_c::Float64  # Critical normal separation (m)
    δ_s_c::Float64  # Critical shear separation (m)
    K_n::Float64    # Normal stiffness (N/m³)
    K_s::Float64    # Shear stiffness (N/m³)
    η::Float64      # Fatigue degradation parameter (per cycle)
end

# Define properties (exaggerated η for fast fatigue)
beam = Material(70e9, 0.3, 0.01, 0.1)
cohesive = CohesiveProperties(50e6, 30e6, 0.001, 0.002, 50e6/0.001, 30e6/0.002, 1e-2)  # η = 1e-2

# Discretization
n_elements = 100
x = range(0, beam.L, length=n_elements)

# Beam deflection and shear displacement (SLB configuration)
function beam_displacements(x, P, beam, a)
    I = beam.h^3 / 12
    w = zeros(length(x))
    u = zeros(length(x))
    
    for i in 1:length(x)
        if x[i] <= a
            w[i] = P * x[i]^2 * (3beam.L - x[i]) / (6 * beam.E * I)
            u[i] = P * x[i]^2 / (2 * beam.E * I) * beam.h / 2
        else
            w[i] = P * beam.L^2 * (3x[i] - beam.L) / (6 * beam.E * I)
            u[i] = P * beam.L^2 / (2 * beam.E * I) * beam.h / 2
        end
    end
    return w, u
end

# Mixed-mode cohesive law (bilinear)
function cohesive_traction(δ_n, δ_s, cohesive, damage)
    σ_max = cohesive.σ_max * (1 - damage)
    τ_max = cohesive.τ_max * (1 - damage)
    K_n = cohesive.K_n * (1 - damage)
    K_s = cohesive.K_s * (1 - damage)
    δ_n_c = cohesive.δ_n_c
    δ_s_c = cohesive.δ_s_c
    
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(δ_n_c^2 + δ_s_c^2)
    
    if δ_eff <= δ_eff_c * σ_max / K_n
        σ = K_n * δ_n
        τ = K_s * δ_s
    elseif δ_eff < δ_eff_c
        λ = (δ_eff_c - δ_eff) / (δ_eff_c - σ_max / K_n)
        σ = σ_max * λ * (δ_n / δ_eff)
        τ = τ_max * λ * (δ_s / δ_eff)
    else
        σ = 0.0
        τ = 0.0
    end
    return σ, τ
end

# Fatigue damage evolution
function update_damage(damage, δ_n, δ_s, cohesive, cycles)
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(cohesive.δ_n_c^2 + cohesive.δ_s_c^2)
    dD_dn = cohesive.η * (δ_eff / δ_eff_c)^2
    new_damage = min(damage + dD_dn * cycles, 1.0)
    return new_damage
end

# SLB simulation with cycle tracking
function simulate_SLB(P, cycles, beam, cohesive)
    a = beam.L / 4  # Initial crack length
    damage = zeros(n_elements)
    
    cycle_points = [0, 200, 400, 600, 800, 1000]
    damage_history = zeros(n_elements, length(cycle_points))
    crack_lengths = zeros(length(cycle_points))
    
    history_idx = 1
    for n in 0:cycles
        w, u = beam_displacements(x, P, beam, a)
        δ_n = w
        δ_s = u
        
        # Update damage and track max values for debugging
        max_δ_eff = 0.0
        for i in 1:n_elements
            δ_eff = sqrt(δ_n[i]^2 + δ_s[i]^2)
            max_δ_eff = max(max_δ_eff, δ_eff)
            if x[i] > a
                damage[i] = update_damage(damage[i], δ_n[i], δ_s[i], cohesive, 1)
            end
        end
        
        # Update crack length
        old_a = a
        for i in 1:n_elements
            if damage[i] >= 1.0 && x[i] > a
                a = x[i]
                break
            end
        end
        
        if n in cycle_points
            damage_history[:, history_idx] .= damage
            crack_lengths[history_idx] = a
            println("Cycle $n: Crack length = $a m, Max δ_eff = $max_δ_eff, Damage change = $(a != old_a)")
            history_idx += 1
        end
    end
    
    return x, damage_history, cycle_points, crack_lengths
end

# Run simulation
P = 10000.0  # Even higher load for significant effect
cycles = 1000
x, damage_hist, cycle_points, crack_lengths = simulate_SLB(P, cycles, beam, cohesive)

# Plot only damage with zoom
xlims_zoom = (0.02, 0.06)  # Focus on crack region
p = plot(xlabel="Position (m)", ylabel="Damage", title="Damage Evolution", xlims=xlims_zoom)
for (i, n) in enumerate(cycle_points)
    plot!(p, x, damage_hist[:, i], label="Cycle $n", linewidth=2)
end
plot(p, size=(800, 400), legend=:topleft)



###########

using Plots
using LinearAlgebra

# Material and geometry properties
struct Material
    E::Float64      # Young's modulus (Pa)
    ν::Float64      # Poisson's ratio
    h::Float64      # Thickness of each layer (m)
    L::Float64      # Length of the beam (m)
end

struct CohesiveProperties
    σ_max::Float64  # Maximum normal strength (Pa)
    τ_max::Float64  # Maximum shear strength (Pa)
    δ_n_c::Float64  # Critical normal separation (m)
    δ_s_c::Float64  # Critical shear separation (m)
    K_n::Float64    # Normal stiffness (N/m³)
    K_s::Float64    # Shear stiffness (N/m³)
    η::Float64      # Fatigue degradation parameter (per cycle)
end

# Define properties
beam = Material(70e9, 0.3, 0.01, 0.1)
cohesive = CohesiveProperties(50e6, 30e6, 0.001, 0.002, 50e6/0.001, 30e6/0.002, 1e-2)

# Discretization
n_elements = 100
x = range(0, beam.L, length=n_elements)

# Beam deflection and shear displacement (SLB configuration)
function beam_displacements(x, P, beam, a)
    I = beam.h^3 / 12
    w = zeros(length(x))
    u = zeros(length(x))
    
    for i in 1:length(x)
        if x[i] <= a
            w[i] = P * x[i]^2 * (3beam.L - x[i]) / (6 * beam.E * I)
            u[i] = P * x[i]^2 / (2 * beam.E * I) * beam.h / 2
        else
            w[i] = P * beam.L^2 * (3x[i] - beam.L) / (6 * beam.E * I)
            u[i] = P * beam.L^2 / (2 * beam.E * I) * beam.h / 2
        end
    end
    return w, u
end

# Mixed-mode cohesive law (bilinear)
function cohesive_traction(δ_n, δ_s, cohesive, damage)
    σ_max = cohesive.σ_max * (1 - damage)
    τ_max = cohesive.τ_max * (1 - damage)
    K_n = cohesive.K_n * (1 - damage)
    K_s = cohesive.K_s * (1 - damage)
    δ_n_c = cohesive.δ_n_c
    δ_s_c = cohesive.δ_s_c
    
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(δ_n_c^2 + δ_s_c^2)
    
    if δ_eff <= δ_eff_c * σ_max / K_n
        σ = K_n * δ_n
        τ = K_s * δ_s
    elseif δ_eff < δ_eff_c
        λ = (δ_eff_c - δ_eff) / (δ_eff_c - σ_max / K_n)
        σ = σ_max * λ * (δ_n / δ_eff)
        τ = τ_max * λ * (δ_s / δ_eff)
    else
        σ = 0.0
        τ = 0.0
    end
    return σ, τ
end

# Fatigue damage evolution
function update_damage(damage, δ_n, δ_s, cohesive, cycles)
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(cohesive.δ_n_c^2 + cohesive.δ_s_c^2)
    dD_dn = cohesive.η * (δ_eff / δ_eff_c)^2
    return min(damage + dD_dn * cycles, 1.0)
end

# SLB simulation with cycle tracking
function simulate_SLB(P, cycles, beam, cohesive)
    a = beam.L / 4  # Initial crack length
    damage = zeros(n_elements)
    
    cycle_points = [0, 200, 400, 600, 800, 1000]
    δ_n_history = zeros(n_elements, length(cycle_points))
    δ_s_history = zeros(n_elements, length(cycle_points))
    σ_history = zeros(n_elements, length(cycle_points))
    τ_history = zeros(n_elements, length(cycle_points))
    damage_history = zeros(n_elements, length(cycle_points))
    crack_lengths = zeros(length(cycle_points))
    
    history_idx = 1
    for n in 0:cycles
        w, u = beam_displacements(x, P, beam, a)
        δ_n = w
        δ_s = u
        
        σ = zeros(n_elements)
        τ = zeros(n_elements)
        for i in 1:n_elements
            σ[i], τ[i] = cohesive_traction(δ_n[i], δ_s[i], cohesive, damage[i])
            if x[i] > a
                damage[i] = update_damage(damage[i], δ_n[i], δ_s[i], cohesive, 1)
            end
        end
        
        old_a = a
        for i in 1:n_elements
            if damage[i] >= 1.0 && x[i] > a
                a = x[i]
                break
            end
        end
        
        if n in cycle_points
            δ_n_history[:, history_idx] .= δ_n
            δ_s_history[:, history_idx] .= δ_s
            σ_history[:, history_idx] .= σ
            τ_history[:, history_idx] .= τ
            damage_history[:, history_idx] .= damage
            crack_lengths[history_idx] = a
            println("Cycle $n: Crack length = $a m")
            history_idx += 1
        end
    end
    
    return x, δ_n_history, δ_s_history, σ_history, τ_history, damage_history, cycle_points, crack_lengths
end

# Run simulation
P = 100000.0  # Load sufficient for visible evolution
cycles = 1000
x, δ_n_hist, δ_s_hist, σ_hist, τ_hist, damage_hist, cycle_points, crack_lengths = simulate_SLB(P, cycles, beam, cohesive)

# Plot all variables with zoom
xlims_zoom = (0., 0.06)
p1 = plot(xlabel="Position (m)", ylabel="Normal Separation (m)", title="Normal Separation", xlims=xlims_zoom)
p2 = plot(xlabel="Position (m)", ylabel="Shear Separation (m)", title="Shear Separation", xlims=xlims_zoom)
p3 = plot(xlabel="Position (m)", ylabel="Normal Traction (Pa)", title="Normal Traction", xlims=xlims_zoom)
p4 = plot(xlabel="Position (m)", ylabel="Shear Traction (Pa)", title="Shear Traction", xlims=xlims_zoom)
p5 = plot(xlabel="Position (m)", ylabel="Damage", title="Damage Evolution", xlims=xlims_zoom)

for (i, n) in enumerate(cycle_points)
    plot!(p1, x, δ_n_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p2, x, δ_s_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p3, x, σ_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p4, x, τ_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p5, x, damage_hist[:, i], label="Cycle $n", linewidth=2)
end

plot(p1, p2, p3, p4, p5, layout=(5,1), size=(800, 1000), legend=:topleft)

# Plot cohesive law (stress vs. separation)
p6 = plot(xlabel="Normal Separation (m)", ylabel="Normal Traction (Pa)", title="Cohesive Law: Normal")
p7 = plot(xlabel="Shear Separation (m)", ylabel="Shear Traction (Pa)", title="Cohesive Law: Shear")

for (i, n) in enumerate(cycle_points)
    # Only plot points where traction is non-zero (cohesive zone)
    mask = (σ_hist[:, i] .> 0) .| (τ_hist[:, i] .> 0)
    scatter!(p6, δ_n_hist[mask, i], σ_hist[mask, i], label="Cycle $n", markersize=3)
    scatter!(p7, δ_s_hist[mask, i], τ_hist[mask, i], label="Cycle $n", markersize=3)
end

# Combine all plots
plot(p1, p2, p3, p4, p5, p6, p7, layout=(7,1), size=(800, 1400), legend=:topleft)



###

using Plots
using LinearAlgebra

# Material and geometry properties
struct Material
    E::Float64      # Young's modulus (Pa)
    ν::Float64      # Poisson's ratio
    h::Float64      # Thickness of each layer (m)
    L::Float64      # Length of the beam (m)
end

struct CohesiveProperties
    σ_max::Float64  # Maximum normal strength (Pa)
    τ_max::Float64  # Maximum shear strength (Pa)
    δ_n_c::Float64  # Critical normal separation (m)
    δ_s_c::Float64  # Critical shear separation (m)
    K_n::Float64    # Normal stiffness (N/m³)
    K_s::Float64    # Shear stiffness (N/m³)
    η::Float64      # Fatigue degradation parameter (per cycle)
end

# Define properties
beam = Material(70e9, 0.3, 0.01, 0.1)
cohesive = CohesiveProperties(50e6, 30e6, 0.001, 0.002, 50e6/0.001, 30e6/0.002, 1e-2)

# Discretization
n_elements = 100
x = range(0, beam.L, length=n_elements)

# Beam deflection and shear displacement (SLB configuration)
function beam_displacements(x, P, beam, a)
    I = beam.h^3 / 12
    w = zeros(length(x))
    u = zeros(length(x))
    
    for i in 1:length(x)
        if x[i] <= a
            w[i] = P * x[i]^2 * (3beam.L - x[i]) / (6 * beam.E * I)
            u[i] = P * x[i]^2 / (2 * beam.E * I) * beam.h / 2
        else
            w[i] = P * beam.L^2 * (3x[i] - beam.L) / (6 * beam.E * I)
            u[i] = P * beam.L^2 / (2 * beam.E * I) * beam.h / 2
        end
    end
    return w, u
end

# Mixed-mode cohesive law (bilinear)
function cohesive_traction(δ_n, δ_s, cohesive, damage)
    σ_max = cohesive.σ_max * (1 - damage)
    τ_max = cohesive.τ_max * (1 - damage)
    K_n = cohesive.K_n * (1 - damage)
    K_s = cohesive.K_s * (1 - damage)
    δ_n_c = cohesive.δ_n_c
    δ_s_c = cohesive.δ_s_c
    
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(δ_n_c^2 + δ_s_c^2)
    
    if δ_eff <= δ_eff_c * σ_max / K_n
        σ = K_n * δ_n
        τ = K_s * δ_s
    elseif δ_eff < δ_eff_c
        λ = (δ_eff_c - δ_eff) / (δ_eff_c - σ_max / K_n)
        σ = σ_max * λ * (δ_n / δ_eff)
        τ = τ_max * λ * (δ_s / δ_eff)
    else
        σ = 0.0
        τ = 0.0
    end
    return σ, τ
end

# Fatigue damage evolution
function update_damage(damage, δ_n, δ_s, cohesive, cycles)
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(cohesive.δ_n_c^2 + cohesive.δ_s_c^2)
    dD_dn = cohesive.η * (δ_eff / δ_eff_c)^2
    return min(damage + dD_dn * cycles, 1.0)
end

# SLB simulation with cycle tracking
function simulate_SLB(P, cycles, beam, cohesive)
    a = beam.L / 4
    damage = zeros(n_elements)
    
    cycle_points = [0, 200, 400, 600, 800, 1000]
    δ_n_history = zeros(n_elements, length(cycle_points))
    δ_s_history = zeros(n_elements, length(cycle_points))
    σ_history = zeros(n_elements, length(cycle_points))
    τ_history = zeros(n_elements, length(cycle_points))
    damage_history = zeros(n_elements, length(cycle_points))
    crack_lengths = zeros(length(cycle_points))
    
    history_idx = 1
    for n in 0:cycles
        w, u = beam_displacements(x, P, beam, a)
        δ_n = w
        δ_s = u
        
        σ = zeros(n_elements)
        τ = zeros(n_elements)
        for i in 1:n_elements
            σ[i], τ[i] = cohesive_traction(δ_n[i], δ_s[i], cohesive, damage[i])
            if x[i] > a
                damage[i] = update_damage(damage[i], δ_n[i], δ_s[i], cohesive, 1)
            end
        end
        
        old_a = a
        for i in 1:n_elements
            if damage[i] >= 1.0 && x[i] > a
                a = x[i]
                break
            end
        end
        
        if n in cycle_points
            δ_n_history[:, history_idx] .= δ_n
            δ_s_history[:, history_idx] .= δ_s
            σ_history[:, history_idx] .= σ
            τ_history[:, history_idx] .= τ
            damage_history[:, history_idx] .= damage
            crack_lengths[history_idx] = a
            println("Cycle $n: Crack length = $a m")
            history_idx += 1
        end
    end
    
    return x, δ_n_history, δ_s_history, σ_history, τ_history, damage_history, cycle_points, crack_lengths
end

# Run simulation
P = 10000.0
cycles = 1000
x, δ_n_hist, δ_s_hist, σ_hist, τ_hist, damage_hist, cycle_points, crack_lengths = simulate_SLB(P, cycles, beam, cohesive)

# Plot spatial evolution
xlims_zoom = (0.02, 0.06)
p1 = plot(xlabel="Position (m)", ylabel="Normal Separation (m)", title="Normal Separation", xlims=xlims_zoom)
p2 = plot(xlabel="Position (m)", ylabel="Shear Separation (m)", title="Shear Separation", xlims=xlims_zoom)
p3 = plot(xlabel="Position (m)", ylabel="Normal Traction (Pa)", title="Normal Traction", xlims=xlims_zoom)
p4 = plot(xlabel="Position (m)", ylabel="Shear Traction (Pa)", title="Shear Traction", xlims=xlims_zoom)
p5 = plot(xlabel="Position (m)", ylabel="Damage", title="Damage Evolution", xlims=xlims_zoom)

for (i, n) in enumerate(cycle_points)
    plot!(p1, x, δ_n_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p2, x, δ_s_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p3, x, σ_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p4, x, τ_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p5, x, damage_hist[:, i], label="Cycle $n", linewidth=2)
end

# Plot cohesive law (stress vs. separation)
p6 = plot(xlabel="Normal Separation (m)", ylabel="Normal Traction (Pa)", title="Cohesive Law: Normal")
p7 = plot(xlabel="Shear Separation (m)", ylabel="Shear Traction (Pa)", title="Cohesive Law: Shear")

for (i, n) in enumerate(cycle_points)
    # Only plot points where traction is non-zero (cohesive zone)
    mask = (σ_hist[:, i] .> 0) .| (τ_hist[:, i] .> 0)
    scatter!(p6, δ_n_hist[mask, i], σ_hist[mask, i], label="Cycle $n", markersize=3)
    scatter!(p7, δ_s_hist[mask, i], τ_hist[mask, i], label="Cycle $n", markersize=3)
end

# Combine all plots
plot(p1, p2, p3, p4, p5, p6, p7, layout=(7,1), size=(800, 1400), legend=:topleft)


###
using Plots
using LinearAlgebra

# Material and geometry properties
struct Material
    E::Float64      # Young's modulus (Pa)
    ν::Float64      # Poisson's ratio
    h::Float64      # Thickness of each layer (m)
    L::Float64      # Length of the beam (m)
end

struct CohesiveProperties
    σ_max::Float64  # Maximum normal strength (Pa)
    τ_max::Float64  # Maximum shear strength (Pa)
    δ_n_c::Float64  # Critical normal separation (m)
    δ_s_c::Float64  # Critical shear separation (m)
    K_n::Float64    # Normal stiffness (N/m³)
    K_s::Float64    # Shear stiffness (N/m³)
    η::Float64      # Fatigue degradation parameter (per cycle)
end

# Define properties
beam = Material(70e9, 0.3, 0.01, 0.1)
cohesive = CohesiveProperties(50e6, 30e6, 0.001, 0.002, 50e6/0.001, 30e6/0.002, 1e-2)

# Discretization
n_elements = 100
x = range(0, beam.L, length=n_elements)

# Beam deflection and shear displacement (SLB configuration)
function beam_displacements(x, P, beam, a)
    I = beam.h^3 / 12
    w = zeros(length(x))
    u = zeros(length(x))
    
    for i in 1:length(x)
        if x[i] <= a
            w[i] = P * x[i]^2 * (3beam.L - x[i]) / (6 * beam.E * I)
            u[i] = P * x[i]^2 / (2 * beam.E * I) * beam.h / 2
        else
            w[i] = P * beam.L^2 * (3x[i] - beam.L) / (6 * beam.E * I)
            u[i] = P * beam.L^2 / (2 * beam.E * I) * beam.h / 2
        end
    end
    return w, u
end

# Mixed-mode cohesive law (bilinear)
function cohesive_traction(δ_n, δ_s, cohesive, damage)
    σ_max = cohesive.σ_max * (1 - damage)
    τ_max = cohesive.τ_max * (1 - damage)
    K_n = cohesive.K_n * (1 - damage)
    K_s = cohesive.K_s * (1 - damage)
    δ_n_c = cohesive.δ_n_c
    δ_s_c = cohesive.δ_s_c
    
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(δ_n_c^2 + δ_s_c^2)
    
    if δ_eff <= δ_eff_c * σ_max / K_n
        σ = K_n * δ_n
        τ = K_s * δ_s
    elseif δ_eff < δ_eff_c
        λ = (δ_eff_c - δ_eff) / (δ_eff_c - σ_max / K_n)
        σ = σ_max * λ * (δ_n / δ_eff)
        τ = τ_max * λ * (δ_s / δ_eff)
    else
        σ = 0.0
        τ = 0.0
    end
    return σ, τ
end

# Fatigue damage evolution
function update_damage(damage, δ_n, δ_s, cohesive, cycles)
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(cohesive.δ_n_c^2 + cohesive.δ_s_c^2)
    dD_dn = cohesive.η * (δ_eff / δ_eff_c)^2
    return min(damage + dD_dn * cycles, 1.0)
end

# SLB simulation with cycle tracking
function simulate_SLB(P, cycles, beam, cohesive)
    a = beam.L / 4
    damage = zeros(n_elements)
    
    cycle_points = [0, 200, 400, 600, 800, 1000]
    δ_n_history = zeros(n_elements, length(cycle_points))
    δ_s_history = zeros(n_elements, length(cycle_points))
    σ_history = zeros(n_elements, length(cycle_points))
    τ_history = zeros(n_elements, length(cycle_points))
    damage_history = zeros(n_elements, length(cycle_points))
    crack_lengths = zeros(length(cycle_points))
    
    history_idx = 1
    for n in 0:cycles
        w, u = beam_displacements(x, P, beam, a)
        δ_n = w
        δ_s = u
        
        σ = zeros(n_elements)
        τ = zeros(n_elements)
        for i in 1:n_elements
            σ[i], τ[i] = cohesive_traction(δ_n[i], δ_s[i], cohesive, damage[i])
            if x[i] > a
                damage[i] = update_damage(damage[i], δ_n[i], δ_s[i], cohesive, 1)
            end
        end
        
        old_a = a
        for i in 1:n_elements
            if damage[i] >= 1.0 && x[i] > a
                a = x[i]
                break
            end
        end
        
        if n in cycle_points
            δ_n_history[:, history_idx] .= δ_n
            δ_s_history[:, history_idx] .= δ_s
            σ_history[:, history_idx] .= σ
            τ_history[:, history_idx] .= τ
            damage_history[:, history_idx] .= damage
            crack_lengths[history_idx] = a
            println("Cycle $n: Crack length = $a m, Max σ = $(maximum(σ)), Max τ = $(maximum(τ))")
            history_idx += 1
        end
    end
    
    return x, δ_n_history, δ_s_history, σ_history, τ_history, damage_history, cycle_points, crack_lengths
end

# Run simulation
P = 10000.0
cycles = 1000
x, δ_n_hist, δ_s_hist, σ_hist, τ_hist, damage_hist, cycle_points, crack_lengths = simulate_SLB(P, cycles, beam, cohesive)

# Plot spatial evolution
xlims_zoom = (0.02, 0.06)
p1 = plot(xlabel="Position (m)", ylabel="Normal Separation (m)", title="Normal Separation", xlims=xlims_zoom)
p2 = plot(xlabel="Position (m)", ylabel="Shear Separation (m)", title="Shear Separation", xlims=xlims_zoom)
p3 = plot(xlabel="Position (m)", ylabel="Normal Traction (Pa)", title="Normal Traction", xlims=xlims_zoom)
p4 = plot(xlabel="Position (m)", ylabel="Shear Traction (Pa)", title="Shear Traction", xlims=xlims_zoom)
p5 = plot(xlabel="Position (m)", ylabel="Damage", title="Damage Evolution", xlims=xlims_zoom)

for (i, n) in enumerate(cycle_points)
    plot!(p1, x, δ_n_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p2, x, δ_s_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p3, x, σ_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p4, x, τ_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p5, x, damage_hist[:, i], label="Cycle $n", linewidth=2)
end

# Plot cohesive law (stress vs. separation)
p6 = plot(xlabel="Normal Separation (m)", ylabel="Normal Traction (Pa)", title="Cohesive Law: Normal")
p7 = plot(xlabel="Shear Separation (m)", ylabel="Shear Traction (Pa)", title="Cohesive Law: Shear")

for (i, n) in enumerate(cycle_points)
    # Relaxed mask: include all points with non-zero separation
    mask = (δ_n_hist[:, i] .> 0) .| (δ_s_hist[:, i] .> 0)
    if any(mask)
        scatter!(p6, δ_n_hist[mask, i], σ_hist[mask, i], label="Cycle $n", markersize=3)
        scatter!(p7, δ_s_hist[mask, i], τ_hist[mask, i], label="Cycle $n", markersize=3)
    else
        println("Cycle $n: No points to plot (all separations zero)")
    end
end

# Set reasonable limits for cohesive law plots
plot!(p6, xlims=(0, 0.0015), ylims=(0, 10e7))
plot!(p7, xlims=(0, 0.0025), ylims=(0, 3e7))

# Combine all plots
plot(p1, p2, p3, p4, p5, p6, p7, layout=(7,1), size=(800, 1400), legend=:topleft)



####

using Plots
using LinearAlgebra

# Material and geometry properties
struct Material
    E::Float64      # Young's modulus (Pa)
    ν::Float64      # Poisson's ratio
    h::Float64      # Thickness of each layer (m)
    L::Float64      # Length of the beam (m)
end

struct CohesiveProperties
    σ_max::Float64  # Maximum normal strength (Pa)
    τ_max::Float64  # Maximum shear strength (Pa)
    δ_n_c::Float64  # Critical normal separation (m)
    δ_s_c::Float64  # Critical shear separation (m)
    K_n::Float64    # Normal stiffness (N/m³)
    K_s::Float64    # Shear stiffness (N/m³)
    η::Float64      # Fatigue degradation parameter (per cycle)
end

# Define properties (adjusted for slower damage)
beam = Material(70e9, 0.3, 0.01, 0.1)
cohesive = CohesiveProperties(50e6, 30e6, 0.002, 0.004, 50e6/0.002, 30e6/0.004, 5e-3)  # Increased δ_c, reduced η

# Discretization
n_elements = 100
x = range(0, beam.L, length=n_elements)

# Beam deflection and shear displacement (SLB configuration)
function beam_displacements(x, P, beam, a)
    I = beam.h^3 / 12
    w = zeros(length(x))
    u = zeros(length(x))
    
    for i in 1:length(x)
        if x[i] <= a
            w[i] = P * x[i]^2 * (3beam.L - x[i]) / (6 * beam.E * I)
            u[i] = P * x[i]^2 / (2 * beam.E * I) * beam.h / 2
        else
            w[i] = P * beam.L^2 * (3x[i] - beam.L) / (6 * beam.E * I)
            u[i] = P * beam.L^2 / (2 * beam.E * I) * beam.h / 2
        end
    end
    return w, u
end

# Mixed-mode cohesive law (bilinear)
function cohesive_traction(δ_n, δ_s, cohesive, damage)
    σ_max = cohesive.σ_max * (1 - damage)
    τ_max = cohesive.τ_max * (1 - damage)
    K_n = cohesive.K_n * (1 - damage)
    K_s = cohesive.K_s * (1 - damage)
    δ_n_c = cohesive.δ_n_c
    δ_s_c = cohesive.δ_s_c
    
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(δ_n_c^2 + δ_s_c^2)
    
    if δ_eff <= δ_eff_c * σ_max / K_n
        σ = K_n * δ_n
        τ = K_s * δ_s
    elseif δ_eff < δ_eff_c
        λ = (δ_eff_c - δ_eff) / (δ_eff_c - σ_max / K_n)
        σ = σ_max * λ * (δ_n / δ_eff)
        τ = τ_max * λ * (δ_s / δ_eff)
    else
        σ = 0.0
        τ = 0.0
    end
    return σ, τ
end

# Fatigue damage evolution
function update_damage(damage, δ_n, δ_s, cohesive, cycles)
    δ_eff = sqrt(δ_n^2 + δ_s^2)
    δ_eff_c = sqrt(cohesive.δ_n_c^2 + cohesive.δ_s_c^2)
    dD_dn = cohesive.η * (δ_eff / δ_eff_c)^2
    return min(damage + dD_dn * cycles, 1.0)
end

# SLB simulation with cycle tracking
function simulate_SLB(P, cycles, beam, cohesive)
    a = beam.L / 4
    damage = zeros(n_elements)
    
    cycle_points = [0, 200, 400, 600, 800, 1000]
    δ_n_history = zeros(n_elements, length(cycle_points))
    δ_s_history = zeros(n_elements, length(cycle_points))
    σ_history = zeros(n_elements, length(cycle_points))
    τ_history = zeros(n_elements, length(cycle_points))
    damage_history = zeros(n_elements, length(cycle_points))
    crack_lengths = zeros(length(cycle_points))
    
    history_idx = 1
    for n in 0:cycles
        w, u = beam_displacements(x, P, beam, a)
        δ_n = w
        δ_s = u
        
        σ = zeros(n_elements)
        τ = zeros(n_elements)
        for i in 1:n_elements
            σ[i], τ[i] = cohesive_traction(δ_n[i], δ_s[i], cohesive, damage[i])
            if x[i] > a
                damage[i] = update_damage(damage[i], δ_n[i], δ_s[i], cohesive, 1)
            end
        end
        
        old_a = a
        for i in 1:n_elements
            if damage[i] >= 1.0 && x[i] > a
                a = x[i]
                break
            end
        end
        
        if n in cycle_points
            δ_n_history[:, history_idx] .= δ_n
            δ_s_history[:, history_idx] .= δ_s
            σ_history[:, history_idx] .= σ
            τ_history[:, history_idx] .= τ
            damage_history[:, history_idx] .= damage
            crack_lengths[history_idx] = a
            println("Cycle $n: Crack length = $a m, Max σ = $(maximum(σ)), Max τ = $(maximum(τ))")
            history_idx += 1
        end
    end
    
    return x, δ_n_history, δ_s_history, σ_history, τ_history, damage_history, cycle_points, crack_lengths
end

# Run simulation
P = 10000.0
cycles = 1000
x, δ_n_hist, δ_s_hist, σ_hist, τ_hist, damage_hist, cycle_points, crack_lengths = simulate_SLB(P, cycles, beam, cohesive)

# Plot spatial evolution
xlims_zoom = (0.02, 0.06)
p1 = plot(xlabel="Position (m)", ylabel="Normal Separation (m)", title="Normal Separation", xlims=xlims_zoom)
p2 = plot(xlabel="Position (m)", ylabel="Shear Separation (m)", title="Shear Separation", xlims=xlims_zoom)
p3 = plot(xlabel="Position (m)", ylabel="Normal Traction (Pa)", title="Normal Traction", xlims=xlims_zoom)
p4 = plot(xlabel="Position (m)", ylabel="Shear Traction (Pa)", title="Shear Traction", xlims=xlims_zoom)
p5 = plot(xlabel="Position (m)", ylabel="Damage", title="Damage Evolution", xlims=xlims_zoom)

for (i, n) in enumerate(cycle_points)
    plot!(p1, x, δ_n_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p2, x, δ_s_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p3, x, σ_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p4, x, τ_hist[:, i], label="Cycle $n", linewidth=2)
    plot!(p5, x, damage_hist[:, i], label="Cycle $n", linewidth=2)
end

# Plot cohesive law for a point just ahead of the crack tip
p6 = plot(xlabel="Normal Separation (m)", ylabel="Normal Traction (Pa)", title="Cohesive Law: Normal")
p7 = plot(xlabel="Shear Separation (m)", ylabel="Shear Traction (Pa)", title="Cohesive Law: Shear")

for (i, n) in enumerate(cycle_points)
    # Find the point just ahead of the crack tip (first point where damage < 1)
    crack_tip = crack_lengths[i]
    idx = findfirst(x -> x > crack_tip && damage_hist[findfirst(isequal(x), x), i] < 1.0, x)
    if !isnothing(idx)
        δ_n_val = δ_n_hist[idx, i]
        δ_s_val = δ_s_hist[idx, i]
        σ_val = σ_hist[idx, i]
        τ_val = τ_hist[idx, i]
        scatter!(p6, [δ_n_val], [σ_val], label="Cycle $n", markersize=5)
        scatter!(p7, [δ_s_val], [τ_val], label="Cycle $n", markersize=5)
        println("Cycle $n: δ_n = $δ_n_val, σ = $σ_val, δ_s = $δ_s_val, τ = $τ_val at x = $(x[idx])")
    else
        println("Cycle $n: No cohesive zone point found")
    end
end

# Set limits to see the full cohesive law
plot!(p6, xlims=(0, 0.0025), ylims=(0, 1e7))
plot!(p7, xlims=(0, 0.005), ylims=(0, 1e7))

# Combine all plots
plot(p1, p2, p3, p4, p5, p6, p7, layout=(7,1), size=(800, 1400), legend=:topleft)