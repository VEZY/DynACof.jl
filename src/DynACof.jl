module DynACof

import CSV.read
using DataFrames
import Dates
import ForwardDiff.derivative # To find Δ (esat slope)


# helpers:
export GDD,is_missing
# ecophysio helpers:
export rH_to_VPD,esat,esat_slope
export physics_constant

# Main functions:
export Meteorology

include("test.jl")
include("helpers.jl")
include("meteo.jl")
include("constants.jl")
include("ecophysio_helpers.jl")

end # module
