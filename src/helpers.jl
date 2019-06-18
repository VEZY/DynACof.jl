"""
    GDD(30.0,27.0,5.0,27.0)

Compute the daily growing degree days (GDD) using the maximum and minimum daily temperature (and optionally
using the average daily temperature, not recommended).

# Arguments
- `Tmax::Float64`: Maximum daily temperature (Celsius degree)
- `Tmin::Float64`: Minimum daily temperature (Celsius degree)
- `MinTT::Float64`: Minimum temperature threshold, also called base temperature (Celsius degree), default to 5.
- `MaxTT::Float64`: Maximum temperature threshold (Celsius degree), optional, default to NULL
- `Round::Bool`: Boolean. Important: round the result to 2 decimal, with default to `TRUE`.
- `Tmean::Float64`: Optional. Average daily temperature (Celsius degree). Only needed if Tmax and Tmin are missing.

Please keep in mind that this function gives an approximation of the degree days. GDD are
usually computed as the integral of hourly (or less) values.
The round argument is provided for convenience, as growing temperatures with less than 3 digits are likely to be
within the measurement error, and will probably have no visible effect on plant phenology.
Caution, use Tmean only if Tmax and Tmin are not available because it tends to give less powerful estimation.


# Return
GDD: Growing degree days (Celsius degree)

# Examples
```julia
# Growing degree days over 10 days :
julia> GDD(30.0,27.0,5.0,27.0)
```

"""
function GDD(Tmax::Float64,Tmin::Float64,MinTT::Float64=5.0,MaxTT::Float64=30.0)::Float64
 Tmean= (Tmax+Tmin)/2.0
 GDD(Tmean,MinTT,MaxTT)
end


function GDD(Tmean::Float64,MinTT::Float64=5.0,MaxTT::Float64=30.0)::Float64
  DD= Tmean-MinTT
  if DD<0.0 || DD>(MaxTT-MinTT)
    DD= 0.0
  end
  DD
end
