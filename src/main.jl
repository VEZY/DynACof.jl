function DynACof(period=NULL, write_it= F,parallel= TRUE,output_f=".RData",input_path="package",
    output_path= Inpath, simulation_name="DynACof",
    file_name= (constants= "constants.jl",site="site.jl",meteo="meteorology.txt",soil="soil.jl",coffee="coffee.jl",tree="tree.jl"))

    # Read in package parameter files (default): 
    if input_path=="package"
        physics= physics_constant()
        # put other parameters here
    else
        physics= read_constants(filepath)
        # Put reading functions for parameters here
    end
    # Parameters= Import_Parameters(path = Inpath, Names= FileName[-grep("Meteo",FileName)])

end

