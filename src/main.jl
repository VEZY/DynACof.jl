"""
    dynacof(;period::Array{String,1}= ["0000-01-01", "0000-01-02"], input_path="package",output_path= input_path, simulation_name="DynACof",
             file_name= (constants= "constants.jl",site="site.jl",meteo="meteorology.txt",soil="soil.jl",coffee="coffee.jl",tree="tree.jl")

# Dynamic Agroforestry Coffee Crop Model

The DynACof process-based model computes plot-scale Net Primary Productivity, carbon allocation, growth, yield, energy, and water balance 
of coffee plantations according to management, while accounting for spatial effects using metamodels from the 3D process-based model
[MAESPA](https: /  / maespa.github.io / ). The model also uses cohorts for the development of the coffee buds and fruits to better represent fruit
carbon demand distribution along the year.

# Arguments
- `period::Array{String,1}`: A vector of two character string as POSIX dates that correspond to the min and max dates for the desired time
period to be returned. The default value ["0000-01-01", "0000-01-02"] makes the function take the min and max values from the meteorology file.
- `input_path::String`: Path to the input parameter list folder. Default to `"package"`, wich makes DynACof use the package default parameter values.
- `output_path::String`: Path pointing to the folder were the results will be written. Default to `"output_path = input_path"`.
- `simulation_name::String`: Character name of the simulation file name when written (`write == true`). Default: `"DynACof"`.
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
|                              | LE_x                   | MJ m-2 d-1          | System / Coffee / Tree / Soil latent heat                                                 |
|                              | H_x                    | MJ m-2 d-1          | System / Coffee / Tree / Soil sensible heat                                               |
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
|                              | WindSpeed_x            | m s-1               | Wind speed at the center of the layer                                               |
|                              | TairCanopy_x           | deg C               | Air tempetature at the center of the layer                                          |
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
|                              | LeafWaterPotential     | Mpa                 | Coffee leaf water potential                                                         |
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


# Examples
```julia
# A simulation with the default parameter files from the package, and an example meteorology file from the `DynACof.jl_inputs` repository:
file= download("https: /  / raw.githubusercontent.com / VEZY / DynACof.jl_inputs / master / meteorology.txt")
Sim, Meteo, Parameters= dynacof(input_path= dirname(file), file_name= (constants= "package",site="package",meteo=basename(file),soil="package",
coffee="package",tree="package"))
rm(file)
```
"""
function dynacof(;period::Array{String,1}= ["0000-01-01", "0000-01-02"], input_path="package",
                 output_path= input_path, simulation_name="DynACof",
                 file_name= (constants= "constants.jl",site="site.jl",meteo="meteorology.txt",soil="soil.jl",
                             coffee="coffee.jl",tree="tree.jl"))


    Parameters= import_parameters(input_path, file_name)
    Meteo= meteorology(normpath(string(input_path,"/",file_name.meteo)), Parameters, period)
    # Setting up the simulation -----------------------------------------------
    # Number of cycles (rotations) to do over the period (given by the Meteo file):
    
    NCycles= ceil(Int64,(maximum(Meteo.year) - minimum(Meteo.year)) / Parameters.AgeCoffeeMax)

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
    age_day_num= vcat(map((x,y) -> collect(x:(1 / y):(x + 1.0 - 1 / y)),age_year,ndaysYear)...)

    Direction= DataFrame(Cycle= cycle_day, Plot_Age= age_day, Plot_Age_num= age_day_num)

    printstyled("Starting a simulation from $(minimum(Meteo.Date)) to $(maximum(Meteo.Date)) over $NCycles plantation cycle(s) \n",
                bold= true, color= :light_green)

    # Potentially make it parallel here
    # for i in 1:NCycles
    #   if i==1
    #     Sim= mainfun(i,Direction,Meteo,Parameters)
    #   else
    #     Sim= vcat(Sim, mainfun(i,Direction,Meteo,Parameters))
    #   end
    # end
  
    Sim= map(x -> mainfun(x,Direction,Meteo,Parameters), 1:NCycles)
    printstyled("Simulation completed successfully \n", bold= true, color= :light_green)

    return vcat(Sim...), Meteo, Parameters
    # Note: to get those elements with a call, just do: `Sim, Meteo, Parameters= dynacof(...)`
end



function mainfun(cy,Direction,Meteo,Parameters)

  # Initializing the table:
  Sim= Direction[Direction.Cycle .== cy,:]
  Met_c= Meteo[Direction.Cycle .== cy,:]
  
  Init_Sim!(Sim,Met_c,Parameters)
  
  # Compute cumulative degree-days based on previous daily DD from semi-hourly data:
  CumulDegreeDays= cumsum(Met_c.DegreeDays)
  
  ## Bud induction window computation ##
  # Bud induction can start only at F_Tffb degree-days after vegetative growth stops.
  # Source: Rodriguez et al. (2011).
  # The following module finds the vegetative growth end day, and add the F_Tffb parameter
  # (Time of first floral buds, in dd), then find the very first flowering of the year
  # and set the vector BudInitPeriod to TRUE between the two dates. So buds will appear
  # between plant F_Tffb parameter and the first flowering day only.

  # Day of vegetative growth end:
  VegetGrowthEndDay= findall(x -> x == Parameters.DVG2, Met_c.DOY)
  # Temporary variables declaration:
  CumsumRelativeToVeget= Array{Float64,2}(undef, length(VegetGrowthEndDay), length(Met_c.Date))
  CumsumRelativeToBudinit= Array{Float64,2}(undef, length(VegetGrowthEndDay), length(Met_c.Date))

  DateBudinit= zeros(Int64, length(VegetGrowthEndDay))
  DateFFlowering= zeros(Int64, length(VegetGrowthEndDay))
  
  for i in 1:length(VegetGrowthEndDay)
    CumsumRelativeToVeget[i,:]= CumulDegreeDays .- CumulDegreeDays[VegetGrowthEndDay[i]-1]
    # Date of first bud initialisation:
    DateBudinit[i]= findlast(CumsumRelativeToVeget[i,:] .< Parameters.F_Tffb)
    CumsumRelativeToBudinit[i,:]= CumulDegreeDays .- CumulDegreeDays[DateBudinit[i]-1]
    # Minimum date of first bud development end (i.e. without dormancy):
    BudDevelEnd= findlast(CumsumRelativeToBudinit[i,:] .< Parameters.F_buds1) - 1
    # Maximum date of first bud development end (i.e. with maximum dormancy):
    MaxDormancy= findlast(CumsumRelativeToBudinit[i,:] .< Parameters.F_buds2) - 1
    # Cumulative rainfall within the period of potential dormancy:
    CumRain= cumsum(Met_c.Rain[BudDevelEnd:MaxDormancy])
    # Effective (real) day of first buds breaking dormancy:
    BudDormancyBreakDay= BudDevelEnd + sum(CumRain .< Parameters.F_rain) - 1
    # Effective day of first flowers:
    DateFFlowering[i]= findlast(CumsumRelativeToBudinit[i,:] .< CumsumRelativeToBudinit[i,BudDormancyBreakDay] .+
                                Parameters.BudInitEnd)
    # Effective dates between which buds can appear
    Sim.BudInitPeriod[DateBudinit[i]:DateFFlowering[i]] .= true
  end

  Sim.BudInitPeriod[CumulDegreeDays .< Parameters.VF_Flowering] .= false

  # Search for the species specific tree function:
  if Parameters.Tree_Species == "No_Shade"
    Treefun! = No_Shade
  else
    Treefun! = Shade_Tree
  end

  
  Sim.ALS= ALS(Elevation= Parameters.Elevation, SlopeAzimut= Parameters.SlopeAzimut, Slope= Parameters.Slope, RowDistance= Parameters.RowDistance,
               Shade= Parameters.Shade, height_coffee= Parameters.Height_Coffee, Fertilization= Parameters.Fertilization,
               ShadeType= Parameters.ShadeType, CoffeePruning= Parameters.CoffeePruning, 
               df_rain= DataFrame(year= Met_c.year, DOY= Met_c.DOY, Rain= Met_c.Rain))

  # Main Loop -----------------------------------------------------------------------------------

  
  p = Progress(length(Sim.LAI),1)

  for i in 1:length(Sim.LAI)
    #   setTxtProgressBar(pb, i)
    next!(p)
    # Shade Tree computation if any
    Treefun!(Sim,Parameters,Met_c,i)
    # Should output at least APAR_Tree, LAI_Tree, T_Tree, Rn_Tree, H_Tree, LE_Tree (sum of transpiration + leaf evap)

    # Coffea computation:
    # CM is in gC m-2soil, so use C content to transform in dry mass
    Sim.LAI[i]= Sim.CM_Leaf[previous_i(i)]  *  Parameters.SLA  /  1000.0  /  Parameters.CC_Leaf
    Sim.LAIplot[i]= Sim.LAIplot[i] + Sim.LAI[i]

    # Light interception ------------------------------------------------------

    Sim.K_Dif[i]= Parameters.k_Dif
    Sim.K_Dir[i]= Parameters.k_Dir

    #APAR coffee
    Sim.PAR_Trans_Tree[i]= Met_c.PAR[i] - Sim.APAR_Tree[i] # PAR above coffee layer
    Sim.APAR_Dif[i]= max(0.0, (Sim.PAR_Trans_Tree[i] * Met_c.FDiff[i]) * (1.0 - exp(-Sim.K_Dif[i] * Sim.LAI[i])))
    APAR_Dir= max(0.0,(Sim.PAR_Trans_Tree[i] * (1.0 - Met_c.FDiff[i])) * (1.0 - exp(-Sim.K_Dir[i] * Sim.LAI[i])))
    # APAR_Dir is not part of Sim because it can be easily computed by Met_c.PARm2d1-Sim.APAR_Dif
    Sim.APAR[i]= APAR_Dir + Sim.APAR_Dif[i]
    Sim.PAR_Trans[i]= Sim.PAR_Trans_Tree[i] - Sim.APAR[i] # PAR above soil layer

    # soil (+canopy evap) water balance ---------------------------------------

    Soilfun!(Sim,Parameters,Met_c,i)

    # Metamodel Coffee leaf water potential
    Sim.LeafWaterPotential[i]= Parameters.LeafWaterPotential(Sim,Met_c,i)


    # Energy balance ----------------------------------------------------------

    # Transpiration Coffee
    Sim.T_Coffee[i]= Parameters.T_Coffee(Sim,Met_c,i)

    # Plot transpiration
    Sim.T_tot[i]= Sim.T_Tree[i] + Sim.T_Coffee[i]

    # Evapo-Transpiration
    Sim.ETR[i]= Sim.T_tot[i] + Sim.E_Soil[i] + Sim.IntercRevapor[i]

    # Latent and Sensible heat fluxes
    Sim.LE_Plot[i]= Sim.ETR[i] * Parameters.λ
    Sim.LE_Coffee[i]= (Sim.T_Coffee[i] + Sim.IntercRevapor[i] * (Sim.LAI[i] / Sim.LAIplot[i])) * Parameters.λ
    Sim.H_Coffee[i]= Parameters.H_Coffee(Sim,Met_c,i)

    # Coffea layer net radiation
    Sim.Rn_Coffee[i]= Sim.H_Coffee[i] + Sim.LE_Coffee[i]

    # Tree LE and Rn (can not compute them in the Tree function because we need IntercRevapor)
    Sim.LE_Tree[i]= (Sim.T_Tree[i] + Sim.IntercRevapor[i] * (Sim.LAI_Tree[i] / Sim.LAIplot[i])) * Parameters.λ
    Sim.Rn_Tree[i]= Sim.H_Tree[i] + Sim.LE_Tree[i]

    # Total plot energy:
    Sim.H_tot[i]= Sim.H_Coffee[i] + Sim.H_Tree[i] + Sim.H_Soil[i]
    Sim.LE_tot[i]= Sim.LE_Coffee[i] + Sim.LE_Tree[i] + Sim.LE_Soil[i]
    Sim.Rn_tot[i]= Sim.Rn_Coffee[i] + Sim.Rn_Tree[i] + Sim.Rn_Soil[i]

    # Tcanopy Coffee : using bulk conductance if no trees, interlayer conductance if trees
    # Source: Van de Griend and Van Boxel 1989.
    if Sim.Height_Tree[i] > Parameters.Height_Coffee

      Sim.TairCanopy[i]= Sim.TairCanopy_Tree[i] + ((Sim.H_Coffee[i] + Sim.H_Soil[i]) * Parameters.MJ_to_W) / 
        (air_density(Sim.TairCanopy_Tree[i],Met_c.Pressure[i] / 10.0) *  Parameters.cp * 
           G_interlay(Wind= Met_c.WindSpeed[i], ZHT = Parameters.ZHT, LAI_top= Sim.LAI_Tree[i], LAI_bot= Sim.LAI[i],
                      Z_top= Sim.Height_Tree[i], extwind = Parameters.extwind))
      Sim.Tleaf_Coffee[i]= Sim.TairCanopy[i] + (Sim.H_Coffee[i] * Parameters.MJ_to_W) / 
        (air_density(Sim.TairCanopy[i],Met_c.Pressure[i] / 10.0) *  Parameters.cp * 
           Gb_h(Wind = Met_c.WindSpeed[i], wleaf= Parameters.wleaf, LAI_lay=Sim.LAI[i], LAI_abv=Sim.LAI_Tree[i],
                ZHT = Parameters.ZHT, Z_top = Sim.Height_Tree[i], extwind= Parameters.extwind))

    else

      Sim.TairCanopy[i]= Met_c.Tair[i] + ((Sim.H_Coffee[i] + Sim.H_Soil[i]) * Parameters.MJ_to_W) / 
        (air_density(Met_c.Tair[i], Met_c.Pressure[i] / 10.0) *  Parameters.cp * 
           G_bulk(Wind = Met_c.WindSpeed[i], ZHT = Parameters.ZHT, Z_top = Parameters.Height_Coffee,
                  LAI = Sim.LAI[i], extwind = Parameters.extwind))

      Sim.Tleaf_Coffee[i]= Sim.TairCanopy[i]+(Sim.H_Coffee[i] * Parameters.MJ_to_W) / 
        (air_density(Sim.TairCanopy[i], Met_c.Pressure[i] / 10.0) *  Parameters.cp  * 
           Gb_h(Wind= Met_c.WindSpeed[i], wleaf= Parameters.wleaf, LAI_lay= Sim.LAI[i], LAI_abv= Sim.LAI_Tree[i],
                ZHT= Parameters.ZHT, Z_top= Parameters.Height_Coffee, extwind= Parameters.extwind))
    end
    # NB: if no trees, TairCanopy_Tree= Tair

    # Recomputing soil temperature knowing TairCanopy

    Sim.TSoil[i]= Sim.TairCanopy[i]+(Sim.H_Soil[i] * Parameters.MJ_to_W) / 
      (air_density(Sim.TairCanopy[i], Met_c.Pressure[i] / 10.0) *  Parameters.cp * 
         G_soilcan(Wind= Met_c.WindSpeed[i], ZHT=Parameters.ZHT, Z_top= max(Sim.Height_Tree[i], Parameters.Height_Coffee),
                   LAI = Sim.LAI_Tree[i] + Sim.LAI[i], extwind= Parameters.extwind))

    Sim.DegreeDays_Tcan[i]= GDD(Sim.Tleaf_Coffee[i], Parameters.MinTT, Parameters.MaxTT)
    
    # Metamodel LUE coffee:
    Sim.lue[i]= Parameters.lue(Sim,Met_c,i)

    #GPP Coffee
    Sim.GPP[i]= Sim.lue[i] * Sim.APAR[i]

    # Maintenance respiration -------------------------------------------------

    # Rm is computed at the beginning of the day on the drymass of the previous day.
    # This is considered as the highest priority for the plant (to maintain its dry mass)

    after_2= i <= 2 ? 0 : 1 # used to start respiration after two days so there are some dry mass.
    # Resprout (branches) wood:
    Sim.Rm_Shoot[i]=  after_2 * (Parameters.pa_Shoot * Sim.DM_Shoot[previous_i(i)] * Parameters.NC_Shoot * 
                      Parameters.MRN *  Parameters.Q10_Shoot^((Sim.TairCanopy[i] - Parameters.TMR) / 10.0))

    # Stump and Coarse roots (perennial wood):
    Sim.Rm_SCR[i]= after_2 *  (Parameters.pa_SCR * Sim.DM_SCR[previous_i(i)] *  Parameters.NC_SCR * Parameters.MRN * 
         Parameters.Q10_SCR^((Sim.TairCanopy[i] - Parameters.TMR) / 10.0))

    # Fruits:
    Sim.Rm_Fruit[i]= after_2 * (Parameters.pa_Fruit * Sim.DM_Fruit[previous_i(i)] *  Parameters.NC_Fruit * Parameters.MRN * 
         Parameters.Q10_Fruit^((Sim.TairCanopy[i] - Parameters.TMR) / 10.0))

    # Leaves:
    Sim.Rm_Leaf[i]= after_2 * (Parameters.pa_Leaf * Sim.DM_Leaf[previous_i(i)] * Parameters.NC_Leaf * Parameters.MRN * 
         Parameters.Q10_Leaf^((Sim.TairCanopy[i] - Parameters.TMR) / 10.0))

    # Fine roots:
    Sim.Rm_FRoot[i]= after_2 * (Parameters.pa_FRoot * Sim.DM_FRoot[previous_i(i)] * Parameters.NC_FRoot * Parameters.MRN * 
         Parameters.Q10_FRoot^((Sim.TairCanopy[i] - Parameters.TMR) / 10.0))

    # Total plant maintenance respiration
    Sim.Rm[i]= Sim.Rm_Fruit[i] + Sim.Rm_Leaf[i] + Sim.Rm_Shoot[i] + Sim.Rm_SCR[i] + Sim.Rm_FRoot[i]


    # Coffee Allocation -------------------------------------------------------

    # Potential use of reserves:
    Sim.Consumption_RE[i]= Parameters.kres * Sim.CM_RE[previous_i(i)]

    # Supply function
    Sim.Supply[i]= max(Sim.GPP[i] - Sim.Rm[i] + Sim.Consumption_RE[i], 0.0)

    # If the respiration is greater than the GPP + reserves use, then take this carbon
    # from mortality of each compartments' biomass equally (not for fruits or reserves):
    Sim.Carbon_Lack_Mortality[i]= -min(0.0, Sim.GPP[i] - Sim.Rm[i] + Sim.Consumption_RE[i])

    # 1-Resprout wood ---------------------------------------------------------
    # Allocation priority 1, see Charbonnier 2012.
    Sim.Alloc_Shoot[i]= Parameters.lambda_Shoot * Sim.Supply[i]
    Sim.NPP_Shoot[i]= Sim.Alloc_Shoot[i] / Parameters.epsilon_Shoot
    Sim.Rg_Shoot[i]= Sim.Alloc_Shoot[i] - Sim.NPP_Shoot[i]
    Sim.Mnat_Shoot[i]= Sim.CM_Shoot[previous_i(i)] / Parameters.lifespan_Shoot
    # Pruning
    if (Sim.Plot_Age[i] >= Parameters.MeanAgePruning) & (Met_c.DOY[i] == Parameters.D_pruning)
      Sim.Mprun_Shoot[i]= Sim.CM_Shoot[previous_i(i)] * Parameters.WoodPruningRate
    end

    Sim.Mortality_Shoot[i]= min((Sim.Mnat_Shoot[i] + Sim.Mprun_Shoot[i]), Sim.CM_Shoot[previous_i(i)])

    # 2-Stump and coarse roots (perennial wood) ------------------------------
    Sim.Alloc_SCR[i]= Parameters.lambda_SCR * Sim.Supply[i]
    Sim.NPP_SCR[i]= Sim.Alloc_SCR[i] / Parameters.epsilon_SCR
    Sim.Rg_SCR[i]= Sim.Alloc_SCR[i] - Sim.NPP_SCR[i]
    Sim.Mnat_SCR[i]= Sim.CM_SCR[previous_i(i)] / Parameters.lifespan_SCR
    Sim.Mortality_SCR[i]= Sim.Mnat_SCR[i]

    # Ratio of number of new nodes per LAI unit as affected by canopy air temperature
    # according to Drinnan & Menzel, 1995
    # Source "0 Effect T on yield and vegetative growth.xlsx", sheet
    # "Std20dComposWinterNodeperBr"
    # NB: computed at the end of the vegetatitve growth only to have Tcan of the
    # whole period already computed
    # NB2 : This is the total number of productive nodes on the coffee plant, i.e. the
    # number of green wood nodes that potentially carry flower buds. Green wood mass (and
    # so number of nodes) are related to leaf area (new leaves appear on nodes) :
    # GUTIERREZ et al. (1998)
    if Met_c.DOY[i] == Parameters.DVG2
      T_VG= Sim.Tleaf_Coffee[(Met_c.year .== Met_c.year[i]) .& (Met_c.DOY .>= Parameters.DVG1) .& (Met_c.DOY .<= Parameters.DVG2)]
      T_VG= sum(T_VG)/length(T_VG)
      Sim.ratioNodestoLAI[Met_c.year .>= Met_c.year[i]] .= Parameters.RNL_base * CN(T_VG)
    end

    # Flower Buds + Flower + Fruits -------------------------------------------

    # (1) Buds induction
    # Buds start appearing for the very first time from 5500 dd. After that,
    # they appear every "Parameters.F_Tffb" degree days until flowering starts
    if Sim.BudInitPeriod[i]
      Sim.Budinit[i]= (Parameters.a_bud+Parameters.b_bud * (Sim.PAR_Trans_Tree[i] / Parameters.FPAR)) *
                       Sim.LAI[i-1] * Sim.ratioNodestoLAI[i-1] * Sim.DegreeDays_Tcan[i]
      # NB: Number of nodes= Sim.LAI[i-1] * Sim.ratioNodestoLAI[i-1]
      Sim.Bud_available[i]= Sim.Budinit[i]
    end
    # NB: number of fruits ~1200  /  year  /  coffee tree, source : Castro-Tanzi et al. (2014)
    # Sim%>%group_by(Plot_Age)%>%summarise(N_Flowers= sum(BudBreak))

    # (2) Cumulative degree days experienced by each bud cohort :
    dd_i= cumsum(Sim.DegreeDays_Tcan[i:-1:previous_i(i,1000)])
    
    # (3) Find the window where buds are under dormancy (find the dormant cohorts)
    # Bud develops during F_buds1 (840) degree days after initiation, so they cannot
    # be dormant less than F_buds1 before i. But they can stay under dormancy until
    # F_buds2 dd maximum, so they cannot be older than F_buds2 dd before i.
    
    OldestDormancy= i - (maximum(findall(dd_i .< Parameters.F_buds2)) - 1)
    YoungestDormancy= i - (maximum(findall(dd_i .< Parameters.F_buds1)) - 1)
    # Idem above (reduce the days computed, F_buds2 is ~300 days and F_buds1 ~80-100 days)

    # (4) Test if the condition of minimum required rain for budbreak is met, and if
    # not, which cohort first met the condition (starting from younger to older cohorts):
    CumRain= cumsum(Met_c.Rain[YoungestDormancy:-1:OldestDormancy])
    # (5) Compute the period were all cohorts have encountered all conditions to break
    # dormancy :
    DormancyBreakPeriod= OldestDormancy:(YoungestDormancy - sum(CumRain .< Parameters.F_rain))

    # (6) Temperature effect on bud phenology
    Sim.Temp_cor_Bud[i]= Parameters.Bud_T_correction()(Sim.Tleaf_Coffee[i])
    
    # (7) Bud dormancy break, Source, Drinnan 1992 and Rodriguez et al., 2011 eq. 13
    Sim.pbreak[i]= 1.0 / (1.0 + exp(Parameters.a_p + Parameters.b_p * Sim.LeafWaterPotential[i]))
    # (8) Compute the number of buds that effectively break dormancy in each cohort:
    Sim.BudBreak_cohort[DormancyBreakPeriod] .=
        map(min, Sim.Bud_available[DormancyBreakPeriod], 
                 Sim.Budinit[DormancyBreakPeriod] .* Sim.pbreak[i] .* Sim.Temp_cor_Bud[DormancyBreakPeriod])
    # NB 1: cannot exceed the number of buds of each cohort
    # NB 2: using Budinit and not Bud_available because pbreak is fitted on total bud cohort

    # (9) Remove buds that did break dormancy from the pool of dormant buds
    Sim.Bud_available[DormancyBreakPeriod]= Sim.Bud_available[DormancyBreakPeriod] .- Sim.BudBreak_cohort[DormancyBreakPeriod]

    # (10) Sum the buds that break dormancy from each cohort to compute the total number
    # of buds that break dormancy on day i :
    Sim.BudBreak[i]= min(sum(Sim.BudBreak_cohort[DormancyBreakPeriod]),Parameters.Max_Bud_Break)
    # Rodriguez et al. state that the maximum number of buds that may break dormancy
    # during each dormancy-terminating episode was set to 12 (see Table 1).

    # Fruits :
    FruitingPeriod= i .- findall(dd_i .< Parameters.F_over) .+ 1
    # NB : Fruits that are older than the FruitingPeriod are overripped

    # Demand from each fruits cohort present on the coffee tree (not overriped),
    # same as Demand_Fruit but keeping each value :
    demand_distribution= logistic_deriv.(dd_i[1:length(FruitingPeriod)], Parameters.u_log, Parameters.s_log) .*
                         [0.0 ; Base.diff(dd_i[1:length(FruitingPeriod)])]
    # NB: we use diff because the values are not evenly distributed (it is not grided, e.g. not 1 by 1 increment)
    demand_distribution[demand_distribution .== Inf] .= 0.0
    Demand_Fruit_Cohort_Period = Sim.BudBreak[FruitingPeriod] .* Parameters.DE_opt .* demand_distribution
    # Total C demand of the fruits :
    Sim.Demand_Fruit[i]= sum(Demand_Fruit_Cohort_Period)
    # C supply to Fruits (i.e. what is left from Supply after removing the consumption
    # by previous compartments and Rm):
    Sim.Supply_Fruit[i]= Sim.Supply[i] - Sim.Alloc_Shoot[i] - Sim.Alloc_SCR[i]

    # Total C allocation to all fruits on day i :
    Sim.Alloc_Fruit[i]= min(Sim.Demand_Fruit[i], Sim.Supply_Fruit[i])
    # Allocation to each cohort, relative to each cohort demand :
    if Sim.Demand_Fruit[i] > 0.0
      Rel_DE= Demand_Fruit_Cohort_Period ./ Sim.Demand_Fruit[i]
    else
      Rel_DE= 0.0
    end
    Sim.Alloc_Fruit_Cohort[FruitingPeriod] .= Sim.Alloc_Fruit[i] .* Rel_DE
    Sim.NPP_Fruit_Cohort[FruitingPeriod] .= Sim.Alloc_Fruit_Cohort[FruitingPeriod] ./ Parameters.epsilon_Fruit
    Sim.CM_Fruit_Cohort[FruitingPeriod] .= Sim.CM_Fruit_Cohort[FruitingPeriod] .+ Sim.NPP_Fruit_Cohort[FruitingPeriod]
    Sim.DM_Fruit_Cohort[FruitingPeriod] .= Sim.CM_Fruit_Cohort[FruitingPeriod] ./ Parameters.CC_Fruit
    # Overriped fruits that fall onto the ground (= to mass of the cohort that overripe) :
    Sim.Overriped_Fruit[i]= Sim.CM_Fruit_Cohort[max(minimum(FruitingPeriod) - 1, 1)]
    # Sim.Overriped_Fruit[i]= Sim.CM_Fruit_Cohort[minimum(FruitingPeriod)-1.0] * Parameters.epsilon_Fruit
    # Duration of the maturation of each cohort born in the ith day (in days):
    Sim.Maturation_duration[FruitingPeriod] .= 1:length(FruitingPeriod)
    # Sucrose content of each cohort:
    Sim.SC[FruitingPeriod] .= Sucrose_cont_perc.(Sim.Maturation_duration[FruitingPeriod], Parameters.S_a, Parameters.S_b,
                                                 Parameters.S_x0, Parameters.S_y0)
    # Sucrose mass of each cohort
    Sim.SM[FruitingPeriod] .= Sim.DM_Fruit_Cohort[FruitingPeriod] .* Sim.SC[FruitingPeriod]
    # Harvest maturity:
    Sim.Harvest_Maturity_Pot[i]= sum(Sim.SM[FruitingPeriod]) / sum(Sim.DM_Fruit_Cohort[FruitingPeriod] .* ((Parameters.S_y0 .+ Parameters.S_a) ./ 100.0))
    # NB : here harvest maturity is computed as the average maturity of the cohorts, because
    # all cohorts present in the Coffea are within the 'FruitingPeriod' window.
    # It could be computed as the percentage of cohorts that are fully mature (Pezzopane
    # et al. 2012 say at 221 days after flowering)
    # Optimal sucrose concentration around 8.8% of the dry mass

    Sim.NPP_Fruit[i]= Sim.Alloc_Fruit[i] / Parameters.epsilon_Fruit
    Sim.Rg_Fruit[i]= Sim.Alloc_Fruit[i] - Sim.NPP_Fruit[i]

    # Harvest. Made one day only for now (TODO: make it a period of harvest)

    if Parameters.harvest == "quantity"
      is_harvest= (Sim.Plot_Age[i] >= Parameters.ageMaturity) & 
                  all(Sim.NPP_Fruit[previous_i.(i,0:10)] .< Sim.Overriped_Fruit[previous_i.(i,0:10)]) &
                  (Sim.CM_Fruit[previous_i(i)] > Parameters.Min_Fruit_CM)
      # Made as soon as the fruit dry mass is decreasing for 10 consecutive days.
      # This condition is met when fruit overriping is more important than fruit NPP
      # for 10 days.
      # This option is the best one when fruit maturation is not well known or when the
      # harvest is made throughout several days or weeks with the assumption that fruits
      # are harvested when mature.
    else
      is_harvest= (Sim.Plot_Age[i]>=Parameters.ageMaturity) &
                  (mean(Sim.Harvest_Maturity_Pot[previous_i.(i,0:9)]) < mean(Sim.Harvest_Maturity_Pot[previous_i.(i,10:19)]))
      # Made as soon as the overall fruit maturation is optimal (all fruits are mature)
    end

    if is_harvest
      # Save the date of harvest:
      Sim.Date_harvest[i]= Met_c.DOY[i]
      Sim.Harvest_Fruit[i]= Sim.CM_Fruit[i-1] + Sim.NPP_Fruit[i] - Sim.Overriped_Fruit[i]
      Sim.Harvest_Maturity[i]= Sim.Harvest_Maturity_Pot[i]
      Sim.CM_Fruit[i-1]= 0.0
      Sim.NPP_Fruit[i]= 0.0
      Sim.Overriped_Fruit[i]= 0.0
      Sim.CM_Fruit_Cohort .= zeros(nrow(Sim))
      # RV: could harvest mature fruits only (To do).
    else
      Sim.Harvest_Fruit[i]= 0.0
    end

    # Leaves ------------------------------------------------------------------

    Sim.Supply_Leaf[i]= Parameters.lambda_Leaf_remain * (Sim.Supply[i] - Sim.Alloc_Fruit[i] - Sim.Alloc_Shoot[i] - Sim.Alloc_SCR[i])

    Sim.Alloc_Leaf[i]= min(Parameters.DELM * (Parameters.Stocking_Coffee / 10000.0) * ((Parameters.LAI_max - Sim.LAI[i]) /
                            (Sim.LAI[i] + Parameters.LAI_max)), 
                           Sim.Supply_Leaf[i])

    Sim.NPP_Leaf[i]= Sim.Alloc_Leaf[i] / Parameters.epsilon_Leaf
    Sim.Rg_Leaf[i]= Sim.Alloc_Leaf[i] - Sim.NPP_Leaf[i]
    Sim.Mnat_Leaf[i]= Sim.CM_Leaf[previous_i(i)] / Parameters.lifespan_Leaf
    Sim.NPP_RE[i]= Sim.NPP_RE[i] + (Sim.Supply_Leaf[i] - Sim.Alloc_Leaf[i])

    Sim.M_ALS[i]= after_2 * max(0.0, Sim.CM_Leaf[previous_i(i)] * Sim.ALS[i])

    if (Sim.Plot_Age[i]>= Parameters.MeanAgePruning) & (Met_c.DOY[i] == Parameters.D_pruning)
      Sim.Mprun_Leaf[i]= Sim.CM_Leaf[previous_i(i)] * Parameters.LeafPruningRate
    else
      Sim.Mprun_Leaf[i]= 0.0
    end

    Sim.Mortality_Leaf[i]= Sim.Mnat_Leaf[i] + Sim.Mprun_Leaf[i] + Sim.M_ALS[i]

    # Fine Roots --------------------------------------------------------------

    Sim.Supply_FRoot[i]= Parameters.lambda_FRoot_remain * (Sim.Supply[i] - Sim.Alloc_Fruit[i] - Sim.Alloc_Shoot[i] - Sim.Alloc_SCR[i])
    Sim.Alloc_FRoot[i]=max(0.0, min(Sim.Alloc_Leaf[i], Sim.Supply_FRoot[i]))
    Sim.NPP_FRoot[i]= Sim.Alloc_FRoot[i] / Parameters.epsilon_FRoot
    Sim.Rg_FRoot[i]= Sim.Alloc_FRoot[i] - Sim.NPP_FRoot[i]
    Sim.NPP_RE[i]= Sim.NPP_RE[i] + (Sim.Supply_FRoot[i] - Sim.Alloc_FRoot[i])
    Sim.Mnat_FRoot[i]= Sim.CM_FRoot[previous_i(i)] / Parameters.lifespan_FRoot
    Sim.Mprun_FRoot[i]= Parameters.m_FRoot * Sim.Mprun_Leaf[i]
    Sim.Mortality_FRoot[i]= Sim.Mnat_FRoot[i] + Sim.Mprun_FRoot[i]

    # Biomass -----------------------------------------------------------------

    CM_tot= Sim.CM_Leaf[previous_i(i)] + Sim.CM_Shoot[previous_i(i)] + Sim.CM_SCR[previous_i(i)] + Sim.CM_FRoot[previous_i(i)]

    Sim.CM_Leaf[i]= Sim.CM_Leaf[previous_i(i)] + Sim.NPP_Leaf[i] - Sim.Mortality_Leaf[i] - 
                    Sim.Carbon_Lack_Mortality[i] * Sim.CM_Leaf[previous_i(i)] / CM_tot
    Sim.CM_Shoot[i]= Sim.CM_Shoot[previous_i(i)] + Sim.NPP_Shoot[i]-Sim.Mortality_Shoot[i] -
                     Sim.Carbon_Lack_Mortality[i] * Sim.CM_Shoot[previous_i(i)] / CM_tot
    Sim.CM_Fruit[i]= Sim.CM_Fruit[previous_i(i)]+ Sim.NPP_Fruit[i] - Sim.Overriped_Fruit[i]
    Sim.CM_SCR[i]= Sim.CM_SCR[previous_i(i)] + Sim.NPP_SCR[i] - Sim.Mortality_SCR[i] -
                   Sim.Carbon_Lack_Mortality[i] * Sim.CM_SCR[previous_i(i)] / CM_tot
    Sim.CM_FRoot[i]= Sim.CM_FRoot[previous_i(i)] + Sim.NPP_FRoot[i] - Sim.Mortality_FRoot[i] -
                     Sim.Carbon_Lack_Mortality[i] * Sim.CM_FRoot[previous_i(i)] / CM_tot
    Sim.CM_RE[i]= Sim.CM_RE[previous_i(i)] + Sim.NPP_RE[i] - Sim.Consumption_RE[i]

    Sim.DM_Leaf[i]= Sim.CM_Leaf[i] / Parameters.CC_Leaf
    Sim.DM_Shoot[i]= Sim.CM_Shoot[i] / Parameters.CC_Shoot
    Sim.DM_Fruit[i]= Sim.CM_Fruit[i] / Parameters.CC_Fruit
    Sim.DM_SCR[i]= Sim.CM_SCR[i] / Parameters.CC_SCR
    Sim.DM_FRoot[i]= Sim.CM_FRoot[i] / Parameters.CC_FRoots
    Sim.DM_RE[i]=Sim.CM_RE[i] / Parameters.CC_SCR

    # Total Respiration and NPP -----------------------------------------------

    Sim.Rg[i]= Sim.Rg_Fruit[i] + Sim.Rg_Leaf[i] + Sim.Rg_Shoot[i]+Sim.Rg_SCR[i] + Sim.Rg_FRoot[i]
    Sim.Ra[i]=Sim.Rm[i] + Sim.Rg[i]
    Sim.NPP[i]= Sim.NPP_Shoot[i] + Sim.NPP_SCR[i] + Sim.NPP_Fruit[i] + Sim.NPP_Leaf[i] + Sim.NPP_FRoot[i]

  end

  return Sim

end