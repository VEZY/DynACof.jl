"""
    dynacof(;period::Array{String,1}= ["0000-01-01", "0000-01-02"], input_path="package",
             file_name= (constants= "constants.jl",site="site.jl",meteo="meteorology.txt",soil="soil.jl",coffee="coffee.jl",tree="tree.jl")

# Dynamic Agroforestry Coffee Crop Model

The DynACof process-based model computes plot-scale Net Primary Productivity, carbon allocation, growth, yield, energy, and water balance
of coffee plantations according to management, while accounting for spatial effects using metamodels from the 3D process-based model
[MAESPA](https://maespa.github.io/). The model also uses cohorts for the development of the coffee buds and fruits to better represent fruit
carbon demand distribution along the year.

# Arguments
- `period::Array{String,1}`: A vector of two character string as POSIX dates that correspond to the min and max dates for the desired time
period to be returned. The default value ["0000-01-01", "0000-01-02"] makes the function take the min and max values from the meteorology file.
- `input_path::String`: Path to the input parameter list folder. Default to `"package"`, wich makes DynACof use the package default parameter values.
- `file_name::NamedTuple{(:constants, :site, :meteo, :soil, :coffee, :tree),NTuple{6,String}}`: A list of input file names :

    + **constants**: Physical constants file. Default: "constants.jl". More info in the corresponding structure: [`constants`](@ref).
    + **site**: Site parameters file name, see details. Default: "site.jl". More info in the corresponding structure: [`site`](@ref)
    + **meteo**: Meteorology file name, see details section. Default: "meteorology.txt". More info in the meteorology reading function [`meteorology`](@ref).
    + **Soil**: Soil parameters file name, see details. Default: "soil.jl". More info in the corresponding structure: [`soil`](@ref).
    + **Coffee**: Coffee parameters file name, see details. Default: "coffee.jl". More info in the corresponding structure: [`coffee`](@ref).
    + **Tree**: Shade tree parameters file name, see details. Default: "tree.jl". More info in the corresponding structure: [`tree`](@ref).

Default input files are provided with the package as an example parameterization. To use the default parameters, you can either set input_path="package" for
using all defaults, or set the desired default file to "package" in `file_name`, *e.g.* to use the default constants, but user-defined other parameters:
`file_name= (constants= "package",site="site.jl",meteo="meteorology.txt",soil="soil.jl",coffee="coffee.jl",tree="tree.jl")`


# Return

Return a three objects Sim, Meteo and Parameters. To get the objects from a dynacof call: `Sim, Meteo, Parameters= dynacof(...)`

- Sim: A data.frame of the simulation outputs at daily time-step:

| Type                         | Var                    | unit                | Definition                                                                          |
|------------------------------|------------------------|---------------------|-------------------------------------------------------------------------------------|
| General                      | Cycle                  | -                   | Plantation cycle ID                                                                 |
|                              | date                   | Posix date (Y-m-d)  | Simulation date                                                                     |
|                              | year                   | Year                | Simulation year                                                                     |
|                              | Plot_Age               | year                | Plantation age (starting at 1)                                                      |
|                              | Plot_Age_num           | year (numeric)      | Numeric age of plantation                                                           |
|                              | LAIplot                | m2 leaves m-2 soil  | Plot (Coffee + Shade Tree if any) Leaf Area Index                                   |
| Suffixes for Coffee organs   | x_RE                   | -                   | Reserves                                                                            |
|                              | x_SCR                  | -                   | Stump and Coarse roots                                                              |
|                              | x_Fruit                | -                   | Fruit                                                                               |
|                              | x_Shoot                | -                   | Resprout wood (= branches)                                                          |
|                              | x_FRoot                | -                   | Fine roots                                                                          |
|                              | x_Leaf                 |                     | Leaves                                                                              |
| Suffixes for Shade Tree org. | x_RE_Tree              | -                   | Reserves                                                                            |
|                              | x_Stem_Tree            | -                   | Stem (= trunk)                                                                      |
|                              | x_Branch_Tree          | -                   | Branches                                                                            |
|                              | x_CoarseRoot_Tree      | -                   | Coarse roots                                                                        |
|                              | x_FRoot_Tree           | -                   | Fine roots                                                                          |
|                              | x_Leaf_Tree            |                     | Leaves                                                                              |
| Energy                       | Rn_tot                 | MJ m-2 d-1          | System net radiation                                                                |
|                              | Rn_Tree                | MJ m-2 d-1          | Shade tree net radiation                                                            |
|                              | Rn_Coffee              | MJ m-2 d-1          | Coffee net radiation                                                                |
|                              | Rn_Soil                | MJ m-2 d-1          | Soil net radiation                                                                  |
|                              | Rn_Soil_SW             | MJ m-2 d-1          | Soil net radiation computed using Shuttleworth & Wallace (1985) for reference       |
|                              | LE_x                   | MJ m-2 d-1          | System / Coffee / Tree / Soil latent heat                                           |
|                              | H_x                    | MJ m-2 d-1          | System / Coffee / Tree / Soil sensible heat                                         |
|                              | Q_Soil                 | MJ m-2 d-1          | Soil heat transport                                                                 |
|                              | Transmittance_Tree     | fraction            | Fraction of light transmitted by the shade trees                                    |
|                              | PAR_Trans_Tree         | MJ m-2 d-1          | Light transmitted by the shade trees canopy                                         |
|                              | PAR_Trans              | MJ m-2 d-1          | Light transmitted by the Coffea canopy                                              |
|                              | K_Dir                  | -                   | Direct light extinction coefficient                                                 |
|                              | K_Dif                  | -                   | Diffuse light extinction coefficient                                                |
|                              | APAR                   | MJ m-2 d-1          | Absorbed PAR by the plant                                                           |
|                              | APAR_Dif               | MJ m-2 d-1          | Absorbed diffuse PAR (Direct is APAR-APAR_Dif)                                      |
|                              | lue                    | gC MJ               | Light use efficiency                                                                |
|                              | Tleaf_Coffee           | deg C               | Coffee canopy temperature computed by DynACof                                       |
|                              | TairCanopy_x           | deg C               | Air tempetature at the center of the layer                                          |
|                              | Gb_h_x                 | m s-1               | Coffee / Tree conductance to heat                                                   |
|                              | Gb_air_canopy          | m s-1               | Bulk (no tree) or canopy layer to canopy layer aerodynamic conductance              |
|                              | air_density_x          | kg m-3              | Air density inside the canopy of the tree or the coffee (see [`air_density`](@ref)) |
|                              | WindSpeed_x            | m s-1               | Wind speed at the center of the layer                                               |
|                              | DegreeDays_Tcan        | deg C               | Growing degree days computed using Coffee Canopy Temperature                        |
| Carbon                       | GPP                    | gC m-2 d-1          | Gross primary productivity                                                          |
|                              | Consumption_RE         | gC m-2 d-1          | Daily reserve consumption                                                           |
|                              | Carbon_Lack_Mortality  | gC m-2 d-1          | Mortality from a higher carbon consumption than Supply                              |
|                              | Rm                     | gC m-2 d-1          | Total Coffee maintenance respiration                                                |
|                              | Rm_x                   | gC m-2 d-1          | Maintenance respiration at organ scale                                              |
|                              | Rg                     | gC m-2 d-1          | Total Coffee growth respiration                                                     |
|                              | Rg_x                   | gC m-2 d-1          | Growth respiration at organ scale                                                   |
|                              | Ra                     | gC m-2 d-1          | Coffee layer autotrophic respiration (=Rm+Rg)                                       |
|                              | Demand_x               | gC m-2 d-1          | C demand at organ scale (fruit, leaf and fine root only)                            |
|                              | Alloc_x                | gC m-2 d-1          | C allocation to organ net of Rm (NPP+Rg)                                            |
|                              | Supply                 | gC m-2 d-1          | C supply at the begining of the day at layer scale (GPP+Reserve consumption-Rm)     |
|                              | Supply_x               | gC m-2 d-1          | C supply to organ, net of Rm                                                        |
|                              | NPP                    | gC m-2 d-1          | Net primary productivity at layer scale                                             |
|                              | NPP_x                  | gC m-2 d-1          | Net primary productivity at organ scale                                             |
|                              | Mnat_x                 | gC m-2 d-1          | Organ natural mortality (= due to lifespan)                                         |
|                              | Mprun_x                | gC m-2 d-1          | Organ mortality due to pruning                                                      |
|                              | M_ALS                  | gC m-2 d-1          | Coffee leaf mortality from American Leaf Spot                                       |
|                              | Mortality_x            | gC m-2 d-1          | Total organ mortality                                                               |
|                              | LAI                    | m2 leaves m-2 soil  | Leaf Area Index                                                                     |
|                              | CM_x                   | gC m-2 d-1          | Organ C mass                                                                        |
|                              | DM_x                   | gDM m-2 d-1         | Organ dry mass                                                                      |
| Fruit development            | BudInitPeriod          | boolean             | Bud initiation period (BIP)                                                         |
|                              | Budinit                | Buds d-1            | Total Number of Buds Initiated per day                                              |
|                              | ratioNodestoLAI        | Nodes LAI-1         | Number of fruiting nodes per LAI unit                                               |
|                              | Temp_cor_Bud           | fraction            | Temperature correction factor for bud development                                   |
|                              | pbreak                 | 0-1                 | Daily probability of bud dormancy break                                             |
|                              | BudBreak               | Buds d-1            | Total number of buds breaking dormancy per day                                      |
|                              | SM                     | g m-2 d-1           | Coffee Fruit Sucrose Mass                                                           |
|                              | SC                     | g Sugar gDM         | Coffee Fruit Sucrose Content                                                        |
|                              | Maturation_duration    | days Fruit cohort-1 | Coffee Fruit Total Maturation Duration for each cohort                              |
|                              | Harvest_Maturity_Pot   | Fraction            | Daily average fruit maturity (0-1)                                                  |
|                              | Date_harvest           | day of year         | date of harvest                                                                     |
|                              | Harvest_Fruit          | gC m-2              | Total fruit carbon mass at harvest                                                  |
|                              | Yield_green            | kg ha-1             | Yield of green coffee bean                                                          |
|                              | Harvest_Maturity       | Fraction            | Average fruit maturity at harvest (0-1)                                             |
|                              | Overriped_Fruit        | gC m-2 d-1          | Overriped fruits that fall onto the ground                                          |
| Water                        | IntercMax              | mm                  | Maximum potential rainfall interception by canopy                                   |
|                              | CanopyHumect           | mm                  | Rainfall interception by canopy                                                     |
|                              | Throughfall            | mm                  | Rainfall not intercepted by the canopy, coming to the soil                          |
|                              | SuperficialRunoff      | mm                  | Water runoff from the superficial soil layer                                        |
|                              | ExcessRunoff           | mm                  | Discharge from the superficial soil layer                                           |
|                              | TotSuperficialRunoff   | mm                  | Sum of discharge+ExcessRunoff                                                       |
|                              | InfilCapa              | mm                  | Superficial water infiltration capacity to first layer of soil                      |
|                              | Infiltration           | mm                  | Superficial water infiltration to first layer of soil                               |
|                              | Drain_[1-3]            | mm                  | Water drainage from soil layer 1, 2 or 3                                            |
|                              | WSurfaceRes            | mm                  | Soil water content from the surface layer                                           |
|                              | W_tot                  | mm                  | Total soil profile water content                                                    |
|                              | W_[1-3]                | mm                  | Soil water content from the layer 1, 2 or 3                                         |
|                              | REW_tot                | -                   | Relative extractable water from the soil                                            |
|                              | REW_[1-3]              | -                   | Relative extractable water from the layer 1, 2 or 3                                 |
|                              | EW_tot                 | mm                  | Extractable water from the soil                                                     |
|                              | EW_[1-3]               | mm                  | Extractable water from the layer 1, 2 or 3                                          |
|                              | SWD                    | mm                  | soil water deficit                                                                  |
|                              | RootWaterExtract_[1-3] | mm                  | Root water extraction for soil layer 1 to 3                                         |
|                              | IntercRevapor          | mm                  | Evaporation by canopy                                                               |
|                              | T_x                    | mm                  | Transpiration at system / Coffee / Tree scale                                           |
|                              | E_Soil                 | mm                  | Soil evaporation                                                                    |
|                              | ETR                    | mm                  | System evapotranspiration                                                           |
|                              | SoilWaterPot           | MPa                 | Soil water potential                                                                |
|                              | PSIL_*                 | Mpa                 | Coffee leaf water potential                                                         |
| Special shade tree variables | LA_Tree                | m2 leaves tree-1    | shade tree leaf area                                                                |
|                              | Crown_H_Tree           | m                   | Crown height                                                                        |
|                              | Trunk_H_Tree           | m                   | Trunk height                                                                        |
|                              | Height_Tree            | m                   | Shade tree total height (used for boundary conductance), set to 0 if no shade trees |
|                              | DBH_Tree               | m                   | Diameter at breast height                                                           |
|                              | LAD_Tree               | m2 m-3              | Shade tree Leaf Area Density                                                        |
|                              | CrownRad_Tree          | m                   | Crown radius                                                                        |
|                              | CrownProj_Tree         | m2 crown tree-1     | Crown projection                                                                    |
|                              | Stocking_Tree          | tree m-2            | Shade tree density                                                                  |
|                              | TimetoThin_Tree        | boolean             | Days on which tree is thinned                                                       |
|                              | MThinning_x_Tree       | gc m-2 d-1          | Mortality due to thining at organ scale                                             |

- Meteo: A data.frame of the input meteorology, potentially coming from the output of [Meteorology()]:

| *Var*           | *unit*      | *Definition*                                 | *If missing*                                                       |
|-----------------|-------------|----------------------------------------------|--------------------------------------------------------------------|
| Date            | POSIXct     | Date in POSIXct format                       | Computed from start date parameter, or set a dummy date if missing |
| year            | year        | Year of the simulation                       | Computed from Date                                                 |
| DOY             | day         | day of the year                              | Computed from Date                                                 |
| Rain            | mm          | Rainfall                                     | Assume no rain                                                     |
| Tair            | Celsius     | Air temperature (above canopy)               | Computed from Tmax and Tmin                                        |
| Tmax            | Celsius     | Maximum air temperature during the day       | Required (error)                                                   |
| Tmin            | Celsius     | Minimum air temperature during the day       | Required (error)                                                   |
| RH              | `%`          | Relative humidity                            | Not used, but prefered over VPD for Rn computation                |
| RAD             | MJ m-2 d-1  | Incident shortwave radiation                 | Computed from PAR                                                  |
| Pressure        | hPa         | Atmospheric pressure                         | Computed from VPD, Tair and Elevation, or alternatively from Tair and Elevation. |
| WindSpeed       | m s-1       | Wind speed                                   | Taken as constant: `Parameters -> WindSpeed`                       |
| CO2             | ppm         | Atmospheric CO2 concentration                | Taken as constant: `Parameters -> CO2`                             |
| DegreeDays      | Celsius     | Growing degree days                          | Computed using [`GDD`](@ref)                                             |
| PAR             | MJ m-2 d-1  | Incident photosynthetically active radiation | Computed from RAD                                                  |
| FDiff           | Fraction    | Diffuse light fraction                       | Computed using [`diffuse_fraction`](@ref) using Spitters et al. (1986) formula  |
| VPD             | hPa         | Vapor pressure deficit                       | Computed from RH                                                   |
| Rn              | MJ m-2 d-1  | Net radiation (will be depreciated)          | Computed using [`Rad_net`](@ref) with RH, or VPD                         |
| DaysWithoutRain | day         | Number of consecutive days with no rainfall  | Computed from Rain                                                 |
| Air_Density     | kg m-3      | Air density of moist air (ρ) above canopy | Computed using [`air_density`](@ref)                       |
| ZEN             | radian      | Solar zenithal angle at noon                 | Computed from Date, Latitude, Longitude and Timezone               |

- Parameters: A list of the input parameters (see [`import_parameters`](@ref), [`constants`](@ref), [`soil`](@ref), [`coffee`](@ref), [`tree`](@ref))

# Details

Almost all variables for coffee exist also for shade trees with the suffix `_Tree` after the name of the variable,
**e.g.**: LAI = coffee LAI, LAI_Tree = shade tree LAI. Special shade tree variables (see return section) are only optional, and there may have more
variables upon parameterization because variables can be added in the parameter file for metamodels_tree or Allometries for example.

# Note

For simulations with custom initialisations (*e.g.* at age > 0), or running a simulation day by day, see [`dynacof_i!`](@ref).

# Examples
```julia
# A simulation with the default parameter files from the package, and an example meteorology file from the `DynACof.jl_inputs` repository:
file= download("https://raw.githubusercontent.com/VEZY/DynACof.jl_inputs/master/meteorology.txt")
Sim, Meteo, Parameters= dynacof(input_path= dirname(file), file_name= (constants= "package",site="package",meteo=basename(file),soil="package",
coffee="package",tree="package"))
rm(file)
```
"""
function dynacof(;period::Array{String,1}= ["0000-01-01", "0000-01-02"], input_path="package",
                 file_name= (constants= "constants.jl",site="site.jl",meteo="meteorology.txt",soil="soil.jl",
                             coffee="coffee.jl",tree="tree.jl"))


    Parameters= import_parameters(input_path, file_name)
    Meteo= meteorology(normpath(string(input_path,"/",file_name.meteo)), Parameters, period)
    # Setting up the simulation -----------------------------------------------
    # Number of cycles (rotations) to do over the period (given by the Meteo file):
    years= unique(Meteo.year)

    NCycles= ceil(Int64,length(years) ./ Parameters.AgeCoffeeMax)

    if NCycles==0
      error("Carefull, minimum allowed simulation length is one year")
    end

    # Setting up the simulation with each plantation rotation (cycle) and plantation age (Plot_Age)
    ndaysYear= zeros(Int64, length(years))
    for i in 1:length(years)
        ndaysYear[i]= nrow(Meteo[Meteo.year .== years[i],:])
    end
    # Variables are re-initialized from one to another cycle so each cycle is independant from the others -> mandatory for
    # parallel processing afterwards

    cycle_year= repeat(1:NCycles, inner= Parameters.AgeCoffeeMax)[1:length(years)]
    cycle_day= vcat(map((x,y) -> repeat([x],y),cycle_year,ndaysYear)...)
    age_year= (0:length(ndaysYear)-1) .% Parameters.AgeCoffeeMax .+ 1
    age_day= vcat(map((x,y) -> repeat([x],inner=y),age_year,ndaysYear)...)
    age_day_num= vcat(map((x,y) -> collect(x:(1 / y):(x + 1.0 - 1 / y)),age_year,ndaysYear)...)

    Direction= DataFrame(Cycle= cycle_day, Plot_Age= age_day, Plot_Age_num= age_day_num)

    printstyled("Starting a simulation from $(minimum(Meteo.Date)) to $(maximum(Meteo.Date)) over $NCycles plantation cycle(s) \n",
                bold= true, color= :light_green)

    # Potentially make it parallel here
    Sim= map(x -> mainfun(x,Direction,Meteo,Parameters), 1:NCycles)
    printstyled("Simulation completed successfully \n", bold= true, color= :light_green)

    return vcat(Sim...), Meteo, Parameters
    # Note: to get those elements with a call, just do: `Sim, Meteo, Parameters= dynacof(...)`
end



function mainfun(cy,Direction,Meteo,Parameters)

  # Initializing the table:
  Sim= Direction[Direction.Cycle .== cy,:]
  Met_c= Meteo[Direction.Cycle .== cy,:]

  initialise!(Sim,Met_c,Parameters)
  if Parameters.Stocking_Coffee > 0.0
    bud_init_period!(Sim,Met_c,Parameters)
    Sim.ALS= ALS(
        Elevation= Parameters.Elevation, SlopeAzimut= Parameters.SlopeAzimut, Slope= Parameters.Slope, RowDistance= Parameters.RowDistance,
        Shade= Parameters.Shade, height_coffee= Parameters.Height_Coffee, Fertilization= Parameters.Fertilization,
        ShadeType= Parameters.ShadeType, CoffeePruning= Parameters.CoffeePruning,
        df_rain= DataFrame(year= Met_c.year, DOY= Met_c.DOY, Rain= Met_c.Rain))
  end

  # Main Loop -----------------------------------------------------------------------------------
  p = Progress(length(Sim.LAI),1)

  for i in 1:length(Sim.LAI)
    next!(p)

    energy_water_models!(Sim,Parameters,Met_c,i) # the soil is in here also
    # Shade Tree computation if any
    if Sim.Stocking_Tree[i] > 0.0
      tree_model!(Sim,Parameters,Met_c,i)
    end

    if Parameters.Stocking_Coffee > 0.0
        coffee_model!(Sim,Parameters,Met_c,i)
    end
  end

  Sim[!,:date] .= Met_c.Date
  Sim[!,:year] .= Met_c.year

  if Parameters.Stocking_Coffee > 0.0
    Sim[!,:Yield_green] .= Sim.Harvest_Fruit ./ 1000.0 .* 10000.0 ./ Parameters.CC_Fruit .* Parameters.FtS
  end

  return Sim
end

"""
    dynacof_i!(i,Sim::DataFrame,Met_c::DataFrame,Parameters)

Using DynACof one iteration after another. Allows to run a DynACof simulation with starting at age > 0 with initializations.

# Arguments
- `i`: Either an integer, or a range giving the day of simulation needed. Match the row index, so `i=366` make a simulation
for the 366th row of Sim and Met.
- `Sim::DataFrame`: The simulation DataFrame (see [`dynacof`](@ref)), initialized using [`dynacof_i_init`](@ref);
- `Met_c::DataFrame`: The meteorology DataFrame (see [`meteorology`](@ref)), initialized using [`dynacof_i_init`](@ref)
- `Parameters`: The parameters for the model (see [`import_parameters`](@ref)), initialized using [`dynacof_i_init`](@ref)

# Examples
```julia

# Making a regular simulation using example data:
file= download("https://raw.githubusercontent.com/VEZY/DynACof.jl_inputs/master/meteorology.txt")

# Initialize the simulation:
Sim, Meteo, Parameters= dynacof_i_init(1:365,input_path= dirname(file), file_name= (constants= "package",site="package",meteo=basename(file),soil="package",coffee="package",tree="package"))
rm(file)

Sim2= copy(Sim)
Meteo2= copy(Meteo)
# Changing the value of Tair in the meteorology for day 366 for S2:
Meteo2.Tair[366]= Meteo2.Tair[366]+10.0

# Make a computation for each:
dynacof_i!(366,Sim,Meteo,Parameters)
dynacof_i!(366,Sim2,Meteo2,Parameters)

# Compare the values of e.g. the maitenance respiration:
Sim.Rm[366]
Sim2.Rm[366]

# To run DynACof for several days, use a range for i:
S= dynacof_i(367:nrow(Meteo),Sim,Meteo,Parameters)
# NB: nrow(Meteo) or nrow(Sim) is the maximum length we can simulate. To increase a simulation,
# initialize it with a wider range for the "Period" argument (see [`dynacof_i_init`](@ref)).
```
"""
function dynacof_i!(i,Sim::DataFrame,Met_c::DataFrame,Parameters)

  if minimum(i) > nrow(Sim)
    error("""Index or range requested ('i') exceeds the range of the simulation. Please provide a maximum
             index/range of $(nrow(Sim)).
             -> If you need a wider range, please initialize a longer simulation using `dynacof_i_init()`.
    """)
  end

  p = Progress(length(i),1)

  for j in collect(i)
    next!(p)
    energy_water_models!(Sim,Parameters,Met_c,j) # the soil is in here also
    # Shade Tree computation if any
    if Sim.Stocking_Tree[j] > 0.0
      tree_model!(Sim,Parameters,Met_c,j)
    end
    # Should output at least APAR_Tree, LAI_Tree, T_Tree, Rn_Tree, H_Tree, LE_Tree (sum of transpiration + leaf evap)
    coffee_model!(Sim,Parameters,Met_c,j)

    Sim.Yield_green[j] = Sim.Harvest_Fruit[j] ./ 1000.0 .* 10000.0 ./ Parameters.CC_Fruit .* Parameters.FtS
  end
end



"""
    dynacof_i_init(i,Sim::DataFrame,Met_c::DataFrame,Parameters)

Initialize a DynACof simulation to be used in dynacof_i.

# Arguments
- `i`: A range giving the days that will be simulated for initialization. Should be >= 365 days.
# Arguments
- `period::Array{String,1}`: A vector of two character string as POSIX dates corresponding to the min and max dates for the whole simulation (used
to pre-allocate the simulation `Data.Frame`). It is *not* the days that will be simulated during initialization, but the whole range possible for simulation afterwards.
The default value ["0000-01-01", "0000-01-02"] makes the function take the min and max values from the meteorology file.
- `input_path::String`: Path to the input parameter list folder. Default to `"package"`, wich makes DynACof use the package default parameter values.
- `file_name::NamedTuple{(:constants, :site, :meteo, :soil, :coffee, :tree),NTuple{6,String}}`: A list of input file names :

    + **constants**: Physical constants file. Default: "constants.jl". More info in the corresponding structure: [`constants`](@ref).
    + **site**: Site parameters file name, see details. Default: "site.jl". More info in the corresponding structure: [`site`](@ref)
    + **meteo**: Meteorology file name, see details section. Default: "meteorology.txt". More info in the meteorology reading function [`meteorology`](@ref).
    + **Soil**: Soil parameters file name, see details. Default: "soil.jl". More info in the corresponding structure: [`soil`](@ref).
    + **Coffee**: Coffee parameters file name, see details. Default: "coffee.jl". More info in the corresponding structure: [`coffee`](@ref).
    + **Tree**: Shade tree parameters file name, see details. Default: "tree.jl". More info in the corresponding structure: [`tree`](@ref).

Default input files are provided with the package as an example parameterization. To use the default parameters, you can either set input_path="package" for
using all defaults, or set the desired default file to "package" in `file_name`, *e.g.* to use the default constants, but user-defined other parameters:
`file_name= (constants= "package",site="site.jl",meteo="meteorology.txt",soil="soil.jl",coffee="coffee.jl",tree="tree.jl")`


# Return

Return three objects: Sim, Meteo and Parameters. To get the objects from the call: `Sim, Meteo, Parameters= dynacof_i_init(...)`. See [`dynacof`](@ref) for more details.

# Examples
```julia

# Making a regular simulation using example data:
file= download("https://raw.githubusercontent.com/VEZY/DynACof.jl_inputs/master/meteorology.txt")
Sim, Meteo, Parameters= dynacof_i_init(1:365,input_path= dirname(file), file_name= (constants= "package",site="package",meteo=basename(file),soil="package",coffee="package",tree="package"))
rm(file)
```
"""
function dynacof_i_init(i;period::Array{String,1}= ["0000-01-01", "0000-01-02"], input_path="package",
                            file_name= (constants= "constants.jl",site="site.jl",meteo="meteorology.txt",soil="soil.jl",
                                         coffee="coffee.jl",tree="tree.jl"))

    if minimum(i)!=1
      error("i must start at 1 for initialization")
    end

    if maximum(i)<365
      error("i must be a range of one year minimum (1:365) for initialization")
    end

    Parameters= import_parameters(input_path, file_name)
    Meteo= meteorology(normpath(string(input_path,"/",file_name.meteo)), Parameters, period)

    # Setting up the simulation -----------------------------------------------
    # Number of cycles (rotations) to do over the period (given by the Meteo file):
    years= unique(Meteo.year)

    NCycles= ceil(Int64,length(years) ./ Parameters.AgeCoffeeMax)

    if NCycles==0
      error("Carefull, minimum allowed simulation length is one year")
    end

    # Setting up the simulation with each plantation rotation (cycle) and plantation age (Plot_Age)
    ndaysYear= zeros(Int64, length(years))
    for i in 1:length(years)
        ndaysYear[i]= nrow(Meteo[Meteo.year .== years[i],:])
    end
    # Variables are re-initialized from one to another cycle so each cycle is independant from the others -> mandatory for
    # parallel processing afterwards

    cycle_year= repeat(1:NCycles, inner= Parameters.AgeCoffeeMax)[1:length(years)]
    cycle_day= vcat(map((x,y) -> repeat([x],y),cycle_year,ndaysYear)...)
    age_year= (0:length(ndaysYear)-1) .% Parameters.AgeCoffeeMax .+ 1
    age_day= vcat(map((x,y) -> repeat([x],inner=y),age_year,ndaysYear)...)
    age_day_num= vcat(map((x,y) -> collect(x:(1 / y):(x + 1.0 - 1 / y)),age_year,ndaysYear)...)

    Direction= DataFrame(Cycle= cycle_day, Plot_Age= age_day, Plot_Age_num= age_day_num)

    printstyled("Starting a simulation from $(minimum(Meteo.Date)) to $(Meteo.Date[maximum(i)]) over $NCycles plantation cycle(s) \n",
                bold= true, color= :light_green)

    # Potentially make it parallel here
    Sim_df= map(x -> mainfun(x,Direction[i,:],Meteo[i,:],Parameters), 1:NCycles)
    Sim_df= vcat(Sim_df...)
    printstyled("Simulation completed successfully \n", bold= true, color= :light_green)

    # Making a simulation DataFrame with the same length as the Meteo, and fill it with the initialization values from Sim_df:
    # Initializing the table:
    Sim= Direction
    Met_c= Meteo

    initialise!(Sim,Met_c,Parameters)
    bud_init_period!(Sim,Met_c,Parameters)

    Sim.ALS= ALS(Elevation= Parameters.Elevation, SlopeAzimut= Parameters.SlopeAzimut, Slope= Parameters.Slope, RowDistance= Parameters.RowDistance,
                 Shade= Parameters.Shade, height_coffee= Parameters.Height_Coffee, Fertilization= Parameters.Fertilization,
                 ShadeType= Parameters.ShadeType, CoffeePruning= Parameters.CoffeePruning,
                 df_rain= DataFrame(year= Met_c.year, DOY= Met_c.DOY, Rain= Met_c.Rain))

    Sim[!,:date] .= Met_c.Date
    Sim[!,:year] .= Met_c.year
    Sim[!,:Yield_green] .= Sim.Harvest_Fruit ./ 1000.0 .* 10000.0 ./ Parameters.CC_Fruit .* Parameters.FtS

    Sim[1:nrow(Sim_df),:]= Sim_df[:,:]

    n_i= min(maximum(i)+1,length(Sim.LAI))
    Sim.LAI[n_i]= Sim.CM_Leaf[maximum(i)]*Parameters.SLA/1000.0/Parameters.CC_Leaf

    if Sim.Stocking_Tree[maximum(i)] > 0.0
        Sim.LAI_Tree[n_i]= Sim.DM_Leaf_Tree[maximum(i)]*(Parameters.SLA_Tree/1000.0)
        Sim.LAIplot[n_i]= Sim.LAIplot[n_i] + Sim.LAI_Tree[n_i]
        Sim.Height_Canopy[n_i]= max(Sim.Height_Tree[maximum(i)], Parameters.Height_Coffee)
    end

    Sim.LAIplot[n_i]= Sim.LAIplot[n_i] + Sim.LAI[n_i]

  return Sim, Meteo, Parameters
end
