module Rini

    #greet() = print("Hello World!")

    using SymPy
    using LinearAlgebra
    using SparseArrays
    

    include("AnalytischesModell.jl")
    include("properties.jl")
    include("FCZM.jl")

    export StaticBeam, properties,  beam_element,cohesive_element,assemble_system,fatigue_degradation,simulate_fatigue

    
end # module Rini
