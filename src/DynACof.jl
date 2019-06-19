module DynACof

import CSV.read

export greet, greet2
export GDD, Meteorology

include("test.jl")
include("helpers.jl")
include("meteo.jl")

greet() = print("Hello World!")
end # module
