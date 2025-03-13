module Rini

    using SymPy
    using LinearAlgebra
    using SparseArrays
    using Statistics
    using Optim
    
    
    include("AnalytischesModell.jl")
    include("properties.jl")
    include("FCZM.jl")
    include("FCZMStructures.jl")
    include("FCZMOptimization.jl")
    
    export StaticBeam, properties,  beam_element,cohesive_element,assemble_system,fatigue_degradation,simulate_fatigue
    export Material, CohesiveProperties, FatigueData, TestSetup
    export FCZMPrediction, objective_function, optimize_FCZM

end # module Rini
