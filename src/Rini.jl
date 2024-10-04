module Rini

    #greet() = print("Hello World!")

    using SymPy

    include("SimpleModel.jl")
    include("properties.jl")
    export StaticBeam, properties
end # module Rini
