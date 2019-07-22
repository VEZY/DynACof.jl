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
export virtual_temp,VPD_to_e,dew_point,paliv_dis
# Meteorology helpers
export Rad_ext,diffuse_fraction,pressure_from_elevation,sun_zenithal_angle
export Rad_net,days_without_rain
# Parameters-related functions
export import_parameters
export constants,site,coffee,soil,tree
export struct_to_tuple
export read_constants,read_site,read_soil,read_coffee,read_tree,read_param_file
export CB,LeafWaterPotential,T_Coffee,H_Coffee,lue,Metamodels_soil
export light_extinction_K_Tree,tree_allometries,metamodels_tree
# Main functions:
export Meteorology

include("test.jl")
include("helpers.jl")
include("meteo.jl")
include("parameters_struct.jl")
include("ecophysio_helpers.jl")
include("meteorology_helpers.jl")
include("import_parameters.jl")

end # module
