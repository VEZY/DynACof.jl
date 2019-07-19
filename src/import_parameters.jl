
function import_parameters(path,Names)
    if path == "package"
        paths= repeat(["package"],length(Names))
        paths_names= keys(Names)
        paths = NamedTuple{paths_names}(paths)       
    else
        paths= map(x -> normpath(string(path,"/",x)),Names)
    end
    if Names.tree==""
        Parameters= (constants= read_constants(paths.constants), site= read_site(paths.site), soil= read_soil(paths.soil),
                     coffee= read_coffee(paths.coffee), tree= ())
    else
        Parameters= (constants= read_constants(paths.constants), site= read_site(paths.site), soil= read_soil(paths.soil),
                     coffee= read_coffee(paths.coffee), tree= read_tree(paths.tree))
    end
    Parameters
end


function read_constants(filepath::String="package")
    if basename(filepath)=="package"
        return physics_constant()
    else
        consts= evalfile(filepath)
        # NB: replace by the following code whenever I make this metaprogramming bit working: 
        # b= map((x,y) -> :($x = $y),collect(keys(consts)),collect(values(consts)))
        # return eval(:(physics_constant($(b...))))
        return physics_constant(consts.Cp,consts.epsi,consts.pressure0,consts.FPAR,consts.g,consts.Rd,consts.Rgas,consts.Kelvin,consts.vonkarman,consts.MJ_to_W,consts.Gsc,consts.σ,consts.H2OMW,consts.W_umol,consts.λ,consts.cl,consts.Dheat)
    end
end

function read_site(filepath::String="package")
    if basename(filepath)=="package"
        return site()
    else
        site_param= evalfile(filepath)
        return site(site_param.Location, site_param.Start_Date, site_param.Latitude, site_param.Longitude, site_param.TimezoneCF,
                    site_param.Elevation, site_param.ZHT, site_param.extwind, site_param.albedo)
    end
end

function read_soil(filepath::String="package")
    if basename(filepath)=="package"
        return soil()
    else
        soil_param= evalfile(filepath)
        return soil(soil_param.TotalDepth, soil_param.Wm1, soil_param.Wm2, soil_param.Wm3, soil_param.Wf1, soil_param.Wf2, 
                    soil_param.Wf3, soil_param.EWMtot, soil_param.IntercSlope, soil_param.WSurfResMax, soil_param.fc, 
                    soil_param.alpha, soil_param.fo, soil_param.kB, soil_param.k_Rn, soil_param.Soil_LE_p, soil_param.PSIE, 
                    soil_param.PoreFrac, soil_param.B, soil_param.RootFraction1, soil_param.RootFraction2, soil_param.RootFraction3, 
                    soil_param.REWc, soil_param.Metamodels_soil)
    end
end

function read_coffee(filepath::String="package")
    if basename(filepath)=="package"
        return coffee()
    else
        cof= evalfile(filepath)
        return coffee(cof.Stocking_Coffee,cof.AgeCoffeeMin,cof.AgeCoffeeMax,cof.SLA,cof.wleaf,cof.DELM,cof.LAI_max,cof.Height_Coffee,
                      cof.D_pruning,cof.MeanAgePruning,cof.LeafPruningRate,cof.WoodPruningRate,cof.k_Dif,cof.k_Dir,cof.kres,cof.DVG1,
                      cof.DVG2,cof.MinTT,cof.MaxTT,cof.RNL_base,cof.VF_Flowering,cof.F_buds1,cof.F_buds2,cof.a_bud,cof.b_bud,cof.F_Tffb,
                      cof.a_p,cof.b_p,cof.F_rain,cof.Max_Bud_Break,cof.ageMaturity,cof.BudInitEnd,cof.F_over,cof.u_log,cof.s_log,cof.S_a,
                      cof.S_b,cof.S_x0,cof.S_y0,cof.Optimum_Berry_DM,cof.kscale_Fruit,cof.harvest,cof.Min_Fruit_CM,cof.FtS,cof.lambda_Shoot,
                      cof.lambda_SCR,cof.lambda_Leaf_remain,cof.lambda_FRoot_remain,cof.lifespan_Leaf,cof.lifespan_Shoot,cof.lifespan_SCR,
                      cof.lifespan_FRoot,cof.m_FRoot,cof.CC_Fruit,cof.CC_Leaf,cof.CC_Shoot,cof.CC_SCR,cof.CC_FRoots,cof.epsilon_Fruit,
                      cof.epsilon_Leaf,cof.epsilon_Shoot,cof.epsilon_SCR,cof.epsilon_FRoot,cof.NC_Fruit,cof.NC_Leaf,cof.NC_Shoot,cof.NC_SCR,
                      cof.NC_FRoot,cof.Q10_Fruit,cof.Q10_Leaf,cof.Q10_Shoot,cof.Q10_SCR,cof.Q10_FRoot,cof.TMR,cof.MRN,cof.pa_Fruit,
                      cof.pa_Leaf,cof.pa_Shoot,cof.pa_SCR,cof.pa_FRoot,cof.DE_opt,cof.Bud_T_correction,cof.SlopeAzimut,cof.Slope,
                      cof.RowDistance,cof.Shade,cof.Fertilization,cof.ShadeType,cof.CoffeePruning,cof.LeafWaterPotential,cof.T_Coffee,
                      cof.H_Coffee,cof.lue)
    end
end


function read_tree(filepath::String="package")
    if basename(filepath)=="package"
        return tree()
    else
        tree_par= evalfile(filepath)
        return tree(tree_par.Tree_Species,tree_par.Species_ID,tree_par.StockingTree_treeha1,tree_par.SLA_Tree,tree_par.wleaf_Tree,
        tree_par.DELM_Tree,tree_par.LAI_max_Tree,tree_par.Leaf_fall_rate_Tree,tree_par.Fall_Period_Tree,tree_par.Thin_Age_Tree,
        tree_par.ThinThresh,tree_par.RateThinning_Tree,tree_par.date_Thin_Tree,tree_par.D_pruning_Tree,tree_par.pruningIntensity_Tree,
        tree_par.m_FRoot_Tree,tree_par.Pruning_Age_Tree,tree_par.lambda_Stem_Tree,tree_par.lambda_Branch_Tree,tree_par.lambda_CR_Tree,
        tree_par.lambda_Leaf_Tree,tree_par.lambda_FRoot_Tree,tree_par.kres_max_Tree,tree_par.Res_max_Tree,tree_par.CC_Leaf_Tree,
        tree_par.CC_wood_Tree,tree_par.epsilon_Branch_Tree,tree_par.epsilon_Stem_Tree,tree_par.epsilon_CR_Tree,tree_par.epsilon_Leaf_Tree,
        tree_par.epsilon_FRoot_Tree,tree_par.epsilon_RE_Tree,tree_par.lifespan_Branch_Tree,tree_par.lifespan_Leaf_Tree,
        tree_par.lifespan_FRoot_Tree,tree_par.lifespan_CR_Tree,tree_par.Kh,tree_par.KhExp,tree_par.Kc,tree_par.KcExp,tree_par.MRN_Tree,
        tree_par.NC_Branch_Tree,tree_par.NC_Stem_Tree,tree_par.NC_CR_Tree,tree_par.NC_Leaf_Tree,tree_par.NC_FRoot_Tree,
        tree_par.Q10Branch_Tree,tree_par.Q10Stem_Tree,tree_par.Q10CR_Tree,tree_par.Q10Leaf_Tree,tree_par.Q10FRoot_Tree,
        tree_par.pa_Branch_Tree,tree_par.pa_Stem_Tree,tree_par.pa_CR_Tree,tree_par.pa_Leaf_Tree,tree_par.pa_FRoot_Tree,tree_par.k,
        tree_par.metamodels_tree,tree_par.Allometries)
    end
end


