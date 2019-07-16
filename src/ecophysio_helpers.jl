"""
    rH_to_VPD(0.5,20,"Allen_1998")

Conversion between vapor pressure (e), vapor pressure deficit (VPD), specific humidity (q), and relative humidity (rH).

# Arguments
- `rH::Float64`: Relative humidity (-)
- `Tair::Float64`: Air temperature (°C)
- `formula::String`: (optional) Formula to be used for the calculation of esat and the slope of esat. One of "Sonntag_1990" (Default),
"Alduchov_1996", or "Allen_1998".

# Returns
VPD, the vapor pressure deficit (kPa)

# Examples
```julia
rH_to_VPD(0.5,20.0,"Allen_1998")
```

# Reference
This function is translated from the R package [bigleaf](https://bitbucket.org/juergenknauer/bigleaf/src/master/).
"""
function rH_to_VPD(rH::Float64, Tair::Float64, formula::String = "Sonntag_1990")::Float64
  if rH > 1 || ismissing(rH)
    error("Relative humidity (rH) has to be between 0 and 1, and no missing values are allowed")
  end
  Esat = esat(Tair,formula)
  Esat - rH * Esat
end


"""
    esat(20,"Allen_1998")

Computes the saturation vapor pressure (Esat)

# Arguments
- `Tair::Float64`: Air temperature (°C)
- `formula::String`: (optional) Formula to be used for the calculation of esat and the slope of esat. One of "Sonntag_1990" (Default),
"Alduchov_1996", or "Allen_1998".

# Returns
Esat, the saturation vapor pressure (kPa)

# Examples
```julia
esat(20.0,"Allen_1998")
```

# Reference
This function is translated from the R package [bigleaf](https://bitbucket.org/juergenknauer/bigleaf/src/master/).

"""
function esat(Tair::Float64, formula::String= "Sonntag_1990")::Float64

  if formula == "Sonntag_1990"
    a = 611.2
    b = 17.62
    c = 243.12
  elseif formula == "Alduchov_1996"
    a = 610.94
    b = 17.625
    c = 243.04
  elseif formula == "Allen_1998"
    a = 610.8
    b = 17.27
    c = 237.3
  else
    error("Wrong formula argument. The formula argument should take values of: " *
    "Sonntag_1990, Alduchov_1996 or Allen_1998")
  end

  a * exp((b * Tair)/(c + Tair))/ 1000.0
end


"""
    esat_slope(20,"Allen_1998")

Computes Δ, the slope of the saturation vapor pressure at given air temperature.

# Arguments
- `Tair::Float64`: Air temperature (°C)
- `formula::String`: (optional) Formula to be used for the calculation of esat and the slope of esat. One of "Sonntag_1990" (Default),
"Alduchov_1996", or "Allen_1998".

# Returns
Δ, the slope of the saturation vapor pressure curve at Tair (``kPa\\ K^{-1}``)

# Examples
```julia
esat_slope(20.0,"Allen_1998")
```
# Reference
This function is translated from the R package [bigleaf](https://bitbucket.org/juergenknauer/bigleaf/src/master/).
"""
function esat_slope(Tair::Float64, formula::String= "Sonntag_1990")::Float64

  if formula == "Sonntag_1990"
    a = 611.2
    b = 17.62
    c = 243.12
  elseif formula == "Alduchov_1996"
    a = 610.94
    b = 17.625
    c = 243.04
  elseif formula == "Allen_1998"
    a = 610.8
    b = 17.27
    c = 237.3
  else
    error("Wrong formula argument. The formula argument should take values of: " *
    "Sonntag_1990, Alduchov_1996 or Allen_1998")
  end

  esat_fun(Tair::Float64)= a * exp((b * Tair)/(c + Tair))
  ForwardDiff.derivative(esat_fun,Tair) / 1000.0
end

"""
    virtual_temp(25.0, 1.5, "Sonntag_1990")
Computes the virtual temperature, *i.e.* the temperature at which dry air would have the same density as moist air at its actual temperature.

# Arguments  
- `Tair::Float64`: Air temperature (°C)
- `pressure::Float64`: Atmospheric pressure (kPa)
- `VPD::Float64`: Vapor pressure deficit (kPa)
- `formula::String`: (optional) Formula to be used for the calculation of esat. One of "Sonntag_1990" (Default),
"Alduchov_1996", or "Allen_1998".
- `C_to_K::Float64`: Celsius degree to Kelvin (*e.g.* 273.15)
- `eps::Float64`: Ratio of the molecular weight of water vapor to dry air 

# Note
`C_to_K` and `eps` can be found using `physics_constant()`

# Returns
T_v, the virtual temperature (°C)

# Examples
```julia
dyn_const= physics_constant()
virtual_temp(20.0,1010.0,1.5)
```

"""
function virtual_temp(Tair::Float64, pressure::Float64, VPD::Float64; formula::String="Sonntag_1990",
   C_to_K::Float64=physics_constant().Kelvin, eps::Float64= physics_constant().eps)::Float64
  e = VPD_to_e(VPD, Tair, formula= "Sonntag_1990")
  Tair = Tair + C_to_K
  Tv = Tair/(1 - (1 - eps) * e/pressure)
  Tv - C_to_K
end


"""
    VPD_to_e(1.5, 25.0, "Sonntag_1990")

Computes the vapor pressure (e) from the vapor pressure deficit (VPD) and the air temperature (Tair)

# Arguments
- `VPD::Float64`: Vapor pressure deficit (kPa)
- `Tair::Float64`: Air temperature (°C)
- `formula::String`: (optional) Formula to be used for the calculation of esat. One of "Sonntag_1990" (Default),
"Alduchov_1996", or "Allen_1998".

# Returns
e, the vapor pressure (kPa)

# Examples
```julia
VPD_to_e(1.5, 25.0, formula= "Sonntag_1990")
```

# Reference
This function is translated from the R package [bigleaf](https://bitbucket.org/juergenknauer/bigleaf/src/master/).

"""
function VPD_to_e(VPD::Float64, Tair::Float64; formula::String="Sonntag_1990")::Float64
  esat(Tair, formula) - VPD
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