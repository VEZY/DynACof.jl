"""
    pressure_from_elevation(1000.0, 25.0, 1.5)

Computes the virtual temperature, *i.e.* the temperature at which dry air would have the same density as moist air at its actual temperature.

# Arguments  
- `Tair::Float64`: Air temperature (°C)
- `pressure::Float64`: Atmospheric pressure (kPa)
- `VPD::Float64`: Vapor pressure deficit (kPa)
- `formula::String`: (optional) Formula to be used for the calculation of esat. One of "Sonntag_1990" (Default),
"Alduchov_1996", or "Allen_1998".
- `C_to_K::Float64`: Celsius degree to Kelvin (*e.g.* 273.15)
- `pressure0::Float64`: reference atmospheric pressure at sea level (kPa)
- `Rd::Float64`: gas constant of dry air (``J\\ kg^{-1}\\ K^{-1}``), source : Foken p. 245
- `g::Float64`: gravitational acceleration (``m\\ s^{-2}``)

# Note
`C_to_K` and `eps` can be found using `physics_constant()`

# Returns
The atmospheric pressure (kPa)

# Examples
```julia
pressure_from_elevation(600.0, 25.0, 1.5)
```

"""
function pressure_from_elevation(elev::Float64, Tair::Float64, VPD::Float64; formula::String="Sonntag_1990",
  C_to_K::Float64= physics_constant().Kelvin, pressure0::Float64= physics_constant().pressure0, 
  Rd::Float64= physics_constant().Rd, g::Float64= physics_constant().g)::Float64

  pressure1= pressure0 / exp(g * elev / (Rd * (Tair + C_to_K)))
  Tv_K= virtual_temp(Tair, pressure1, VPD, formula = formula) + C_to_K
  pressure0 / exp(g * elev / (Rd * Tv_K))
end


"""
    diffuse_fraction(DOY::Int64, RAD::Float64, Latitude::Float64; formula::String="Spitters",Gsc::Float64=physics_constant().Gsc)

Computes the daily diffuse fraction from the total daily incident radiation

# Arguments  
- `DOY::Int64`: Day Of Year from 1st January (day)
- `RAD::Float64`: Incident total radiation (MJ m-2 d-1)
- `Latitude::Float64`: Latitude (deg)
- `formula::String`: (Optionnal) Model type, one of `Spitters`, `Page` or `Gopinathan`
- `Gsc::Float64`: (Optionnal) The solar constant (W m-2), default to `physics_constant().Gsc` (= 1367).

# Details 
The daily extra-terrestrial radiation at a plane parallel to the earth surface (`S0` or `H0` depending on the source) is computed following
Khorasanizadeh and Mohammadi (2016).
The daily diffuse fraction is computed following DB models from :

* Spitters et al. (1986): used in de Bilt in Netherlands, stated that their model is 
valid for a wide range of climate conditions  
* Page (1967) using the data from 10 widely-spread sites in the 40N to 40S latitude belt  
* Gopinathan and Soler (1995) from 40 widely distributed locations in the latitude range of 36S to 60N.  

# Note
`C_to_K` and `eps` can be found using `physics_constant()`

# Returns
``Hd/H``: the daily diffuse fraction of light (%)

# References
* Duffie, J.A. and W.A. Beckman, Solar engineering of thermal processes. 2013: John Wiley & Sons.
Gopinathan, K. and A. Soler, Diffuse radiation models and monthly-average, daily, diffuse data for
a wide latitude range. Energy, 1995. 20(7): p. 657-667.  
* Kalogirou, S.A., Solar energy engineering: processes and systems. 2013: Academic Press.
Khorasanizadeh, H. and K. Mohammadi, Diffuse solar radiation on a horizontal surface:
Reviewing and categorizing the empirical models. Renewable and Sustainable Energy Reviews,
2016. 53: p. 338-362.  
* Liu, B.Y.H. and R.C. Jordan, The interrelationship and characteristic distribution of direct,
diffuse and total solar radiation. Solar Energy, 1960. 4(3): p. 1-19.  
* Page, J. The estimation of monthly mean values of daily total short wave radiation on vertical
and inclined surfaces from sunshine records 40S-40N. in Proceedings of the United Nations
Conference on New Sources of Energy: Solar Energy, Wind Power and Geothermal Energy, Rome, Italy. 1967.  
* Spitters, C.J.T., H.A.J.M. Toussaint, and J. Goudriaan, Separating the diffuse and direct
component of global radiation and its implications for modeling canopy photosynthesis Part I.
Components of incoming radiation. Agricultural and Forest Meteorology, 1986. 38(1): p. 217-229.  

# Examples
```julia
# Daily diffuse fraction of january 1st at latitude 35 N, with a RAD of 25 MJ m-2 day-1 :
diffuse_fraction(1,25.0,35.0)
```

"""
function diffuse_fraction(DOY::Int64, RAD::Float64, Latitude::Float64; formula::String="Spitters",
    Gsc::Float64=physics_constant().Gsc)::Float64

  S0= Rad_ext(DOY,Latitude,Gsc)

  if S0<=0
    TRANS= 0.0
  else
    TRANS = RAD/S0
  end 

  if formula=="Spitters"
    if TRANS < 0.07
        FDIF= 1.0
    elseif 0.07 <= TRANS < 0.35
        FDIF= 1.0 - 2.3 * (TRANS-0.07)^2.0
    elseif 0.35 <= TRANS <0.75
        FDIF= 1.33 - 1.46 * TRANS
    else
        FDIF= 0.23
    end
    return FDIF
  elseif formula=="Page"
    return 1.0 - 1.13 * TRANS
  elseif formula=="Gopinathan"
    return 0.91138 - 0.96225 * TRANS
  else
    error("Wrong value for formula argument. It should be one of Spitters, Page or Gopinathan.")
  end
end


"""
    Rad_ext(1000.0, 25.0, 1.5)

Computes the virtual temperature, *i.e.* the temperature at which dry air would have the same density as moist air at its actual temperature.

# Arguments  
- `DOY::Int64`: Ordinal date (integer): day of year from 1st January (day)
- `Latitude::Float64`: Latitude (deg)
- `Gsc::Float64`: The solar constant (W m-2), default to `physics_constant().Gsc` (= 1367).

# Returns
`S0`, the daily extra-terrestrial radiation (``MJ\\ m^{-2}\\ d^{-1}``)


# References

Khorasanizadeh, H. and K. Mohammadi, Diffuse solar radiation on a horizontal surface:
Reviewing and categorizing the empirical models. Renewable and Sustainable Energy Reviews,
2016. 53: p. 338-362.

# Examples
```julia
# Daily extra-terrestrial radiation on january 1st at latitude 35 N :
Rad_ext(1,35.0)
```

"""
function Rad_ext(DOY::Int64,Latitude::Float64,Gsc::Float64=physics_constant().Gsc)::Float64
    solar_declin= 23.45*sin°(((float(DOY)+284.0)*360.0)/365.0)
    sunset_hour_angle= acos°(-tan°(Latitude) * tan°(solar_declin))
    S0= (86400.0/π) * Gsc * (1.0 + 0.033 * cos°(float(360*DOY)/365.0)) * (cos°(Latitude) * 
        cos°(solar_declin) * sin°(sunset_hour_angle) + (π * sunset_hour_angle / 180.0) *
        sin°(Latitude) * sin°(solar_declin))
    S0*10^-6
end


"""
    sun_zenithal_angle(DOY::Int64, Latitude::Float64)

Computes the sun zenithal angle at noon (solar time).

# Arguments  
- `DOY::Int64`: Ordinal date (integer): day of year from 1st January (day)
- `Latitude::Float64`: Latitude (deg)

# Returns
`ZEN`, the sun zenithal angle (`radian`)

# References

[solartime](https://cran.r-project.org/web/packages/solartime/) R package from Thomas Wutzler, and more specificly the
`computeSunPositionDoyHour` function (considering the hour at noon).

# Examples
```julia
# Daily extra-terrestrial radiation on january 1st at latitude 35 N :
sun_zenithal_angle(1,35.0)
```

"""
function sun_zenithal_angle(DOY::Int64, Latitude::Float64)::Float64
  fracYearInRad= 2.0 * π * (DOY - 1.0)/365.24

  SolDeclRad=
  ((0.33281 - 22.984 * cos(fracYearInRad) - 
  0.3499 * cos(2.0 * fracYearInRad) - 0.1398 * cos(3.0 * fracYearInRad) + 
  3.7872 * sin(fracYearInRad) + 0.03205 * sin(2.0 * fracYearInRad) + 
  0.07187 * sin(3.0 * fracYearInRad))/180.0 * π)

SolElevRad= asin(sin(SolDeclRad) * sin(Latitude/180.0 * π) + cos(SolDeclRad) * cos(Latitude/180.0 * π))
acos(sin(SolElevRad))
end