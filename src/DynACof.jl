module DynACof

using DataFrames
using CSV
import Dates
import ForwardDiff.derivative # To find Δ (esat slope)
import Optim.optimize, Optim.minimizer


# Helper functions:
export is_missing
export cos°,sin°,tan°,acos°,asin°,atan°
# ecophysio helpers:
export rH_to_VPD,esat,esat_slope,GDD
export virtual_temp,VPD_to_e,dew_point
# Constants
export physics_constant
# Meteorology helpers
export Rad_ext,diffuse_fraction,pressure_from_elevation,sun_zenithal_angle
export Rad_net,days_without_rain
# Main functions:
export Meteorology

include("test.jl")
include("helpers.jl")
include("meteo.jl")
include("constants.jl")
include("ecophysio_helpers.jl")
include("meteorology_helpers.jl")

end # module
