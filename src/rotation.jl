"""
    rotation(Meteo, rotation_length)

Computes a DataFrame with three columns:

- Cycle: The rotation index
- Plot_Age: The age of the rotation (days)
- Plot_Age_num: The age in a numerical form (0 -> 1)

# Arguments

- `Meteo`: The daily meteo file, *e.g.* output from [`meteorology`](@ref)
- `rotation_length`: the length of the rotation (year)
"""
function rotation(Meteo, rotation_length)
    years = unique(Meteo.year)

    NCycles = ceil(Int64, length(years) ./ rotation_length)

    if NCycles == 0
        error("Carefull, minimum allowed simulation length is one year")
    end

    # Setting up the simulation with each plantation rotation (cycle) and plantation age (Plot_Age)
    ndaysYear = zeros(Int64, length(years))
    for i in 1:length(years)
        ndaysYear[i] = nrow(Meteo[Meteo.year.==years[i], :])
    end
    # Variables are re-initialized from one to another cycle so each cycle is independant from the others -> mandatory for
    # parallel processing afterwards

    cycle_year = repeat(1:NCycles, inner=rotation_length)[1:length(years)]
    cycle_day = vcat(map((x, y) -> repeat([x], y), cycle_year, ndaysYear)...)
    age_year = (0:length(ndaysYear)-1) .% rotation_length .+ 1
    age_day = vcat(map((x, y) -> repeat([x], inner=y), age_year, ndaysYear)...)
    age_day_num = vcat(map((x, y) -> collect(x:(1/y):(x+1.0-1/y)), age_year, ndaysYear)...)

    DataFrame(Cycle=cycle_day, Plot_Age=age_day, Plot_Age_num=age_day_num)
end
