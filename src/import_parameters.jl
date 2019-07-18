
function Import_Parameters(path,Names)

    physics_constant()
    # Parameters= c(,site(),coffee(),soil(),if(!is.null(Names$Tree)){Tree()}else{list(Tree_Species= "No_Shade")})

end


function read_constants(filepath::String="package")
    if filepath=="package"
        return physics_constant()
    else
        consts= evalfile(filepath)
        return physics_constant(consts.Cp,consts.eps,consts.pressure0,consts.FPAR,consts.g,consts.Rd,consts.Rgas,
        consts.Kelvin,consts.vonkarman,consts.MJ_to_W,consts.Gsc,consts.σ,consts.H2OMW,consts.W_umol,consts.λ,consts.cl,consts.Dheat)
    end
end

function read_site()

end
