function dynacof(;period::Array{String,1}= ["0000-01-01", "0000-01-02"], input_path="package",
                 output_path= input_path, simulation_name="DynACof",
                 file_name= (constants= "constants.jl",site="site.jl",meteo="meteorology.txt",soil="soil.jl",
                             coffee="coffee.jl",tree="tree.jl"))


    Parameters= import_parameters(input_path, file_name)
    Meteo= meteorology(normpath(string(input_path,"/",file_name.meteo)), Parameters, period)
    # Setting up the simulation -----------------------------------------------
    # Number of cycles (rotations) to do over the period (given by the Meteo file):
    
    # NCycles= ceiling((max(Meteo$year)-min(Meteo$year))/Parameters$AgeCoffeeMax)
    
    # if NCycles==0
    #   error("Carefull, minimum allowed simulation length is one year")
    # end

end

