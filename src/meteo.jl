"""
Meteorology(file = NULL, period = NULL, Parameters = Import_Parameters())

Import the meteorology data, check its format, and eventually compute missing variables.
# Arguments
- `file::String`: The meteorology file path.
- `period::Array{String,1}`: A vector of two character string as POSIX dates that correspond to the min and max dates for the desired time period to be returned.
The default value ["0000-01-01", "0000-01-02"] makes the function take the min and max values from the meteorology file.
- `Parameters::Date`: A list of parameters
* Start_Date: optional, the Posixct date of the first meteo file record. Only needed if the Date column is missing.
* FPAR      : Fraction of global radiation corresponding to PAR radiation, only needed if either RAD or PAR is missing.
* Elevation : elevation of the site (m), only needed if atmospheric pressure is missing
* Latitude  : latitude of the site (degree), only needed if the diffuse fraction of light is missing
* WindSpeed : constant wind speed (m s-1), only needed if windspeed is missing
* CO2       : constant atmospheric ``CO_2`` concentration (ppm), only needed if ``CO_2`` is missing
* MinTT     : minimum temperature threshold for degree days computing (Celsius), see [GDD()]
* MaxTT     : maximum temperature threshold for degree days computing (Celsius), see [GDD()]
* albedo    : site shortwave surface albedo, only needed if net radiation is missing, see [Rad_net()]

Details: The imported file is expected to be at daily time-step. The albedo is used to compute the system net radiation that is then
used to compute the soil net radiation using an extinction coefficient with the plot LAI following the Shuttleworth & Wallace (1985)
formulation. This computation is likely to be depreciated in the near future as the computation has been replaced by a metamodel. It
is kept for information for the moment.

| *Var*           | *unit*      | *Definition*                                 | *If missing*                                                       |
|-----------------|-------------|----------------------------------------------|--------------------------------------------------------------------|
| Date            | POSIXct     | Date in POSIXct format                       | Computed from start date parameter, or set a dummy date if missing |
| year            | year        | Year of the simulation                       | Computed from Date                                                 |
| DOY             | day         | day of the year                              | Computed from Date                                                 |
| Rain            | mm          | Rainfall                                     | Assume no rain                                                     |
| Tair            | Celsius     | Air temperature (above canopy)               | Computed from Tmax and Tmin                                        |
| Tmax            | Celsius     | Maximum air temperature during the day       | Required (error)                                                   |
| Tmin            | Celsius     | Minimum air temperature during the day       | Required (error)                                                   |
| RH              | `%`          | Relative humidity                            | Not used, but prefered over VPD for Rn computation                 |
| RAD             | MJ m-2 d-1  | Incident shortwave radiation                 | Computed from PAR                                                  |
| Pressure        | hPa         | Atmospheric pressure                         | Computed from VPD, Tair and Elevation, or alternatively from Tair and Elevation. |
| WindSpeed       | m s-1       | Wind speed                                   | Taken as constant: `Parameters -> WindSpeed`                          |
| CO2             | ppm         | Atmospheric CO2 concentration                | Taken as constant: `Parameters -> CO2`                                |
| DegreeDays      | Celsius     | Growing degree days                          | Computed using [GDD()]                                             |
| PAR             | MJ m-2 d-1  | Incident photosynthetically active radiation | Computed from RAD                                                  |
| FDiff           | Fraction    | Diffuse light fraction                       | Computed using [Diffuse_d()] using Spitters et al. (1986) formula  |
| VPD             | hPa         | Vapor pressure deficit                       | Computed from RH                                                   |
| Rn              | MJ m-2 d-1  | Net radiation (will be depreciated)          | Computed using [Rad_net()] with RH, or VPD                         |
| DaysWithoutRain | day         | Number of consecutive days with no rainfall  | Computed from Rain                                                 |
| Air_Density     | kg m-3      | Air density of moist air (``\\rho``) above canopy | Computed using [bigleaf::air.density()]                      |
| ZEN             | radian      | Solar zenithal angle at noon                 | Computed from Date, Latitude, Longitude and Timezone               |

Note: It is highly recommended to set the system environment timezone to the one from the meteorology file. If not, the function try to use the Timezone
from the parameter files to set it. When in doubt, set it to UTC (`Sys.setenv(TZ="UTC")`), as for [Aquiares()].

# Returns
A daily meteorology data.frame (invisibly).


# Examples
```julia
Met_c= Meteorology()
```
"""
function Meteorology(file::String, Parameters::Dict, period::Array{String,1}= ["0000-01-01", "0000-01-02"])
    period_date= Dates.Date.(period_date)

    MetData= CSV.read(file; copycols=true);

    if is_missing(MetData,"Date")
        if !is_missing(Parameters,"Start_Date")
            MetData[:Date] =
            collect(Dates.Date(Parameters["Start_Date"]):Dates.Day(1):
            (Dates.Date(Parameters["Start_Date"]) + Dates.Day(nrow(MetData)-1)))
            warn_var("Date","Start_Date from Parameters","warn")
        end
    else
        MetData[:Date] = collect(Dates.Date("2000-01-01"):Dates.Day(1):
        (Dates.Date(Dates.Date("2000-01-01")) + Dates.Day(nrow(MetData)-1)))
        warn_var("Date","dummy 2000-01-01", "warn")
    end

    if period_date != ["0000-01-01", "0000-01-02"]
        if period_date[1]<min(MetData.Date)| period_date[2]>max(MetData.Date)
            error("Given period is not covered by the meteorology file")
        else
            MetData= MetData[period_date[1] .<= MetData.Date .<= period_date[2], :]
        end
    end

    if is_missing(MetData,"RAD")
        if !is_missing(MetData,"PAR")
            MetData[:RAD] = MetData[:PAR] ./ Parameters["FPAR"]
            warn_var("RAD", "PAR", "warn")
        else
            warn_var("RAD", "PAR", "error")
        end
    end

    if is_missing(MetData,"PAR")
        if !is_missing(MetData,"RAD")
            MetData[:PAR] = MetData[:RAD] .* Parameters["FPAR"]
            warn_var("PAR", "RAD", "warn")
        else
            warn_var("PAR", "RAD", "error")
        end
    end
    MetData.PAR[MetData.PAR.<0.1, :] .= 0.1

    if is_missing(MetData,"Tmin") || is_missing(MetData,"Tmax")
        warn_var("Tmin and/or Tmax","error")
    end

    if is_missing(MetData,"Tair")
        MetData[:Tair] = (MetData.Tmax .+ MetData.Tmin) ./ 2.0
        warn_var("Tair","the equation (MetData.Tmax-MetData.Tmin)/2","warn")
    end

    if is_missing(MetData,"VPD")
        if !is_missing(MetData,"RH")
            MetData[:VPD] = rH_to_VPD.(MetData.RH ./ 100.0, MetData.Tair) .* 10.0 # hPa
            warn_var("VPD","RH and Tair using bigleaf::rH.to.VPD","warn")
        else
            warn_var("VPD","RH","error")
        end
    end

    1+2

    # # Missing air pressure:
    # if(is.null(MetData$Pressure)){
    #   if(!is.null(Parameters$Elevation)){
    #     if(!is.null(MetData$VPD)){
    #       bigleaf::pressure.from.elevation(elev = Parameters$Elevation,
    #                                        Tair = MetData$Tair,
    #                                        VPD = MetData$VPD)*10
    #       # Return in kPa
    #       warn_var("Pressure",
    #                paste("Elevation, Tair and VPD",
    #                                  "using bigleaf::pressure.from.elevation"),
    #                "warn")
    #     }else{
    #       bigleaf::pressure.from.elevation(elev = Parameters$Elevation,
    #                                        Tair = MetData$Tair)*10
    #       # Return in kPa
    #       warn_var("Pressure",
    #                paste("Elevation and Tair",
    #                                  "using bigleaf::pressure.from.elevation"),
    #                "warn")
    #     }
    #   }else{
    #     warn_var("Pressure","Elevation","error")
    #   }
    # }
    #
    # # Missing rain:
    # if(is.null(MetData$Rain)){
    #   MetData$Rain= 0 # assume no rain
    #   warn_var("Rain","constant (= 0, assuming no rain)","warn")
    # }
    #
    # # Missing wind speed:
    # if(is.null(MetData$WindSpeed)){
    #   if(!is.null(Parameters$WindSpeed)){
    #     MetData$WindSpeed= Parameters$WindSpeed # assume constant windspeed
    #     warn_var("WindSpeed","constant (= Parameters$WindSpeed)","warn")
    #   }else{
    #     warn_var("WindSpeed",  "Parameters$WindSpeed (constant value)","error")
    #   }
    # }
    # MetData$WindSpeed[MetData$WindSpeed<0.01]= 0.01
    # # Missing atmospheric CO2 concentration:
    # if(is.null(MetData$CO2)){
    #   if(!is.null(Parameters$CO2)){
    #     MetData$CO2= Parameters$CO2 # assume constant windspeed
    #     warn_var("CO2","constant (= Parameters$CO2)","warn")
    #   }else{
    #     warn_var("CO2",  "Parameters$CO2 (constant value)","error")
    #   }
    # }
    #
    # # Missing DegreeDays:
    # if(is.null(MetData$DegreeDays)){
    #   MetData$DegreeDays=
    #     GDD(Tmax= MetData$Tmax,Tmin= MetData$Tmin, MinTT= Parameters$MinTT,
    #         MaxTT = Parameters$MaxTT)
    #   warn_var("DegreeDays","Tmax, Tmin and MinTT","warn")
    # }
    #
    # # Missing diffuse fraction:
    # if(is.null(MetData$FDiff)){
    #   MetData$FDiff=
    #     Diffuse_d(DOY = MetData$DOY, RAD = MetData$RAD,
    #               Latitude = Parameters$Latitude,type = "Spitters")
    #   warn_var("FDiff","DOY, RAD and Latitude using Diffuse_d()","warn")
    # }
    #
    #
    # MetData$year= lubridate::year(MetData$Date)
    # MetData$DOY= lubridate::yday(MetData$Date)
    #
    # # Correct the noon hour by the Timezone if the user use TZ="UTC":
    # if(Sys.timezone()=="UTC"|Sys.timezone()=="GMT"){
    #   cor_tz= Parameters$TimezoneCF*60*60
    # }else{
    #   # Else R use the user time zone (with warning).
    #   warning("Meteo file uses this time-zone: ",Sys.timezone(),". Set it to "UTC" if you want to use ",
    #           "the timezone from your parameter file")
    #   cor_tz= 1
    # }
    #
    # # Solar zenithal angle at noon (radian):
    # MetData$ZEN=
    #   solartime::computeSunPosition(timestamp = MetData$Date+60*60*12+cor_tz,
    #                                 latDeg = Parameters$Latitude,
    #                                 longDeg = Parameters$Longitude)%>%
    #   as.data.frame()%>%{sin(.$elevation)}%>%acos(.)
    #
    # # Compute net radiation using the Allen et al. (1998) equation :
    #
    # if(!is.null(MetData$RH)){
    #   MetData$Rn= Rad_net(DOY = MetData$DOY,RAD = MetData$RAD,Tmax = MetData$Tmax,
    #                       Tmin = MetData$Tmin, Rh =  MetData$RH,
    #                       Elevation = Parameters$Elevation,Latitude = Parameters$Latitude,
    #                       albedo = Parameters$albedo)
    # }else if(!is.null(MetData$VPD)){
    #   MetData$Rn= Rad_net(DOY = MetData$DOY,RAD = MetData$RAD,Tmax = MetData$Tmax,
    #                       Tmin = MetData$Tmin, VPD =  MetData$VPD,
    #                       Elevation = Parameters$Elevation,Latitude = Parameters$Latitude,
    #                       albedo = Parameters$albedo)
    # }
    #
    # DaysWithoutRain= Rain= NULL # To avoid notes by check
    # MetData= as.data.table(MetData)
    # MetData[, DaysWithoutRain := 0]; MetData[Rain > 0, DaysWithoutRain := 1]
    # MetData$DaysWithoutRain= sequence(MetData[,.N,cumsum(DaysWithoutRain)]$N)-1
    # MetData= as.data.frame(MetData)
    #
    # MetData$Air_Density= bigleaf::air.density(MetData$Tair,MetData$Pressure/10)
    #
    # # Force to keep only the input variable the model need to avoid any issues:
    # Varnames= c('year','DOY','Date','Rain','Tair','RH','RAD','Pressure',
    #             'WindSpeed','CO2','DegreeDays','PAR','FDiff',
    #             'VPD','Rn','Tmax','Tmin','DaysWithoutRain','Air_Density','ZEN')
    # MetData= MetData[colnames(MetData)%in%Varnames]
    # MetData[,-c(1:3)]= round(MetData[,-c(1:3)],4)
    #
    # attr(MetData,"unit")=
    #   data.frame(Varnames,
    #              unit=c("year","day","POSIXct date","mm","Celsius","%","MJ m-2 d-1","hPa",
    #                     "m s-1","ppm","Celsius","MJ m-2 d-1","Fraction","hPa","MJ m-2 d-1",
    #                     "Celsius","Celsius","day","kg m-3","rad"))
    #
    # message("Meteo computation done")
    # message(paste("n", crayon::green$bold$underline("Meteo computation done")))
    return MetData
end
