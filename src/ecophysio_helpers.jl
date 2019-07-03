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
function rH_to_VPD(rH::Float64, Tair::Float64, formula::String = "Sonntag_1990")
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
function esat(Tair::Float64, formula::String= "Sonntag_1990")

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
function esat_slope(Tair::Float64, formula::String)

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
