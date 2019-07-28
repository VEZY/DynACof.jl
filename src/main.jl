function dynacof(;period::Array{String,1}= ["0000-01-01", "0000-01-02"], input_path="package",
                 output_path= input_path, simulation_name="DynACof",
                 file_name= (constants= "constants.jl",site="site.jl",meteo="meteorology.txt",soil="soil.jl",
                             coffee="coffee.jl",tree="tree.jl"))


    Parameters= import_parameters(input_path, file_name)
    Meteo= meteorology(normpath(string(input_path,"/",file_name.meteo)), Parameters, period)
    # Setting up the simulation -----------------------------------------------
    # Number of cycles (rotations) to do over the period (given by the Meteo file):
    
    NCycles= ceil(Int64,(maximum(Meteo.year)-minimum(Meteo.year))/Parameters.AgeCoffeeMax)

    if NCycles==0
      error("Carefull, minimum allowed simulation length is one year")
    end

    # Setting up the simulation with each plantation rotation (cycle) and plantation age (Plot_Age) 
    years= unique(Meteo.year)
    ndaysYear= zeros(Int64, length(years))
    for i in 1:length(years)
        ndaysYear[i]= nrow(Meteo[Meteo.year .== years[i],:])
    end
    # Variables are re-initialized from one to another cycle so each cycle is independant from the others -> mandatory for
    # parallel processing afterwards

    cycle_year= repeat(1:NCycles, inner= Parameters.AgeCoffeeMax)[1:length(unique(Meteo.year))]
    cycle_day= vcat(map((x,y) -> repeat([x],y),cycle_year,ndaysYear)...)
    age_year= collect(Parameters.AgeCoffeeMin:min(Parameters.AgeCoffeeMax,length(years)))
    age_day= vcat(map((x,y) -> repeat([x],inner=y),age_year,ndaysYear)...)
    
    Direction= DataFrame(Cycle= cycle_day, Plot_Age= age_day)

    printstyled("Starting a simulation from $(minimum(Meteo.Date)) to $(maximum(Meteo.Date)) over $NCycles plantation cycle(s) \n",
                bold= true, color= :light_green)
end

