module DynACof

import CSV.read
import DataFrames.DataFrame

export greet, greet2

# helpers:
export GDD, Meteorology,is_missing

include("test.jl")
include("helpers.jl")
include("meteo.jl")

greet() = print("Hello World!")
end # module
