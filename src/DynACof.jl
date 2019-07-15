module DynACof

import CSV.read
import DataFrames
import Dates
import ForwardDiff.derivative # To find Δ (esat slope)


# helpers:
export GDD,is_missing
# ecophysio helpers:
export rH_to_VPD,esat,esat_slope
export constant

# Main functions:
export Meteorology

include("test.jl")
include("helpers.jl")
include("meteo.jl")
include("ecophysio_helpers.jl")

end # module
