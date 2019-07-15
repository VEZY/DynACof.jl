"""
Constants used in the DynACof package

This function defines the following constants:
- Cp: specific heat of air for constant pressure (``J\\ K^{-1}\\ kg^{-1}``), Source: Allen 1998 FAO Eq. 8 p. 32
- pressure0: reference atmospheric pressure at sea level (Pa)
- FPAR: Fraction of global radiation that is PAR (source: MAESPA model)
- g: gravitational acceleration (``m\\ s^{-2}``)
- Rd: gas constant of dry air (``J\\ kg^{-1}\\ K^{-1}``), source : Foken p. 245
- Rgas: universal gas constant (``J\\ mol^{-1}\\ K^{-1}``)
- Kelvin: conversion degree Celsius to Kelvin
- vonkarman: von Karman constant (-)
- MJ_to_W: coefficient to convert MJ into W (``W\\ MJ^{-1}``)
- Gsc: solar constant (``W\\ m^{-2}=J\\ m^{-2}\\ s^{-1}``), source : Khorasanizadeh and Mohammadi (2016)
- σ (sigma): Stefan-Boltzmann constant (``W\\ m^{-2}\\ K^{-4}``)
- H2OMW: Conversion factor from kg to mol for H2O (``kg\\ mol^{-1}``)
- W_umol: Conversion factor from watt to micromole for H2O (``W\\ \\mu mol^{-1}``)
- λ (lambda): Latent heat of vaporization (``MJ\\ kg_{H2O}^{-1}``)
- cl: Drag coefficient per unit leaf area (``m\\ s^{-1}``)
- Dheat: Molecular diffusivity for heat (``m\\ s^{-1}``)

Values are partly burrowed from [bigleaf::bigleaf.constants()](https://www.rdocumentation.org/packages/bigleaf/versions/0.7.0/topics/bigleaf.constants)

# References
- Allen, R. G., et al. (1998). "Crop evapotranspiration-Guidelines for computing crop water requirements-FAO Irrigation and drainage paper 56."  300(9): D05109.
- Foken, T, 2008: Micrometeorology. Springer, Berlin, Germany.
- Khorasanizadeh, H. and K. Mohammadi (2016). "Diffuse solar radiation on a horizontal surface: Reviewing and categorizing the empirical models." Renewable and Sustainable Energy Reviews 53: 338-362.
"""
function physics_constant(;
    cp::Float64= 1013*10^-6,
    pressure0::Float64 = 101325.0,
    FPAR::Float64      = 0.5,
    g::Float64         = 9.81,
    Rd::Float64        = 287.0586,
    Rgas::Float64      = 8.314,
    Kelvin::Float64    = 273.15,
    vonkarman::Float64 = 0.41,
    MJ_to_W::Float64   = 10^-6,
    Gsc::Float64       = 1367.0,            # also found 1366 in Kalogirou (2013)
    σ::Float64         = 5.670367e-08,
    H2OMW::Float64     = 18.e-3,
    W_umol::Float64    = 4.57,
    λ::Float64         = 2.45,
    cl::Float64        = 0.4,
    Dheat::Float64     = 21.5e-6
    )

    Constants(cp,pressure0,FPAR,g,Rd,Rgas,Kelvin,vonkarman,MJ_to_W,Gsc,σ,H2OMW,W_umol,λ,cl,Dheat)
end

struct Constants
    cp::Float64
    pressure0::Float64
    FPAR::Float64
    g::Float64
    Rd::Float64
    Rgas::Float64
    Kelvin::Float64
    vonkarman::Float64
    MJ_to_W::Float64
    Gsc::Float64
    σ::Float64
    H2OMW::Float64
    W_umol::Float64
    λ::Float64
    cl::Float64
    Dheat::Float64
end
