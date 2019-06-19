"""
    GDD(30.0,27.0,5.0,27.0)

Compute the daily growing degree days (GDD) using the maximum and minimum daily temperature.

# Arguments
- `Tmax::Float64`: Maximum daily temperature (Celsius degree)
- `Tmin::Float64`: Minimum daily temperature (Celsius degree)
- `MinTT::Float64`: Minimum temperature threshold, also called base temperature (Celsius degree), default to 5.
- `MaxTT::Float64`: Maximum temperature threshold (Celsius degree), optional, default to 30.0

Please keep in mind that this function gives an approximation of the degree days.
GDD are normally computed as the integral of hourly (or less) values.

# Return
GDD: Growing degree days (Celsius degree)

# Examples
```julia
GDD(30.0,27.0,5.0,27.0)
0.0
```
"""
function GDD(Tmax::Float64,Tmin::Float64,MinTT::Float64=5.0,MaxTT::Float64=30.0)::Float64
 Tmean= (Tmax+Tmin)/2.0
 GDD(Tmean,MinTT,MaxTT)
end

"""
    GDD(25.,5.0,28.0)

Compute the daily growing degree days (GDD) directly from the daily mean
temperature.

# Arguments
- `Tmean::Float64`: Optional. Average daily temperature (Celsius degree).
- `MinTT::Float64`: Minimum temperature threshold, also called base temperature (Celsius degree), default to 5.
- `MaxTT::Float64`: Maximum temperature threshold (Celsius degree), optional, default to 30.0

# Return
GDD: Growing degree days (Celsius degree)

# Examples
```julia
GDD(25.0,5.0,28.0)
20.0
GDD(5.0,5.0,28.0)
0.0
```
"""
function GDD(Tmean::Float64,MinTT::Float64=5.0,MaxTT::Float64=30.0)::Float64
  DD= Tmean-MinTT
  if DD<0.0 || DD>(MaxTT-MinTT)
    DD= 0.0
  end
  DD
end
