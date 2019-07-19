
function Import_Parameters(path,Names)

    read_constants()
    # Parameters= c(,site(),coffee(),soil(),if(!is.null(Names$Tree)){Tree()}else{list(Tree_Species= "No_Shade")})

end


function read_constants(filepath::String="package")
    if filepath=="package"
        return physics_constant()
    else
        consts= evalfile(filepath)
        #  return physics_constant(consts.Cp,consts.epsi,consts.pressure0,consts.FPAR,consts.g,consts.Rd,consts.Rgas,consts.Kelvin,consts.vonkarman,consts.MJ_to_W,consts.Gsc,consts.σ,consts.H2OMW,consts.W_umol,consts.λ,consts.cl,consts.Dheat)
        return consts
    end
end

function read_site(filepath::String="package")
    if filepath=="package"
        return site()
    else
        site_param= evalfile(filepath)
        site_values= site(site_param.Location, site_param.Start_Date, site_param.Latitude, site_param.Longitude, site_param.TimezoneCF,
         site_param.Elevation, site_param.ZHT, site_param.extwind, site_param.albedo)
        return site_values
    end
end
