
struct Material
    E::Float64          # Young's modulus (Pa)
    ν::Float64          #   Poisson's ratio
    thickness::Float64  # Thickness (m)
end

# CZM Eigenschaften
struct CohesiveProperties
    K0::Float64     # Initial stiffness (N/m³)
    σ_max::Float64  # Maximum traction (Pa)
    δ_max::Float64  # Maximum separation (m)
    δ_c::Float64    # Critical separation (m)
    m::Float64      # Fatigue degradation parameter
    G_c:: Float64   # Critical energy release rate (J/m2)

    
    function CohesiveProperties(σ_max::Float64,δ_max::Float64, δ_c::Float64, m::Float64)
        K0  = σ_max / δ_max
        G_c = σ_max* δ_c/2
        return new(K0, σ_max,δ_max, δ_c, m,G_c)
    end
end

struct FatigueData
    crack_length::Float64
    separation::Float64
end

struct TestSetup
    L_steel::Float64  
    L_layered::Float64 
    width_layered::Float64
    width_steel::Float64  
    n_elem::Int 
    dx::Float64  
    nodes::Vector{Float64}
    start_pos::Float64  
    end_pos::Float64  
    start_elem::Int
    end_elem::Int
    n_elem_layered::Int
end
