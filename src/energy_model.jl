"""
Coffee crop model

Computes the coffee crop growth and yield. This function is called from [`dynacof`](@ref) and should not be called
by the user.

# Return

Nothing, modify the DataFrame of simulation `Sim` in place. See [`dynacof`](@ref) for more details.

"""
function energy_water_models!(Sim,Parameters,Met_c,i)
    # NB: carefull the order of execution is important here, some function need previous ones in particular order
    if Sim.Stocking_Tree[i] > 0.0
        light_model_tree!(Sim,Parameters,Met_c,i)
    end
    
    light_model_coffee!(Sim,Parameters,Met_c,i)

    if Sim.Stocking_Tree[i] > 0.0
        energy_model_tree!(Sim,Parameters,Met_c,i)
    end
    soil_model!(Sim,Parameters,Met_c,i)
    energy_model_coffee!(Sim,Parameters,Met_c,i)
    
    # Soil temperature: we have to know TairCanopy to compute it, but we have to know H_Soil in energy_model_coffee!
    # so we have to compute the soil before the coffee. 
    Sim.TSoil[i]= Sim.TairCanopy[i] + (Sim.H_Soil[i] * Parameters.MJ_to_W) / 
                  (air_density(Sim.TairCanopy[i], Met_c.Pressure[i]/10.0) * Parameters.cp *
                   G_soilcan(Wind= Met_c.WindSpeed[i], ZHT=Parameters.ZHT, Z_top= max(Sim.Height_Tree[i], Parameters.Height_Coffee),
                             LAI = Sim.LAI_Tree[i]  +  Sim.LAI[i], extwind= Parameters.extwind))

    balance_model!(Sim,Parameters,Met_c,i) # Energy balance
end


function light_model_tree!(Sim,Parameters,Met_c,i)
    Base.invokelatest(Parameters.k, Sim,Met_c,i)
  
    Sim.APAR_Dif_Tree[i]= (Met_c.PAR[i] * Met_c.FDiff[i]) * (1.0-exp(-Sim.K_Dif_Tree[i] * Sim.LAI_Tree[i]))
    Sim.APAR_Dir_Tree[i]= (Met_c.PAR[i]*(1-Met_c.FDiff[i])) * (1.0-exp(-Sim.K_Dir_Tree[i]*Sim.LAI_Tree[i]))
  
    Sim.APAR_Tree[i]= max(0.0,Sim.APAR_Dir_Tree[i]+Sim.APAR_Dif_Tree[i])
  
    Sim.Transmittance_Tree[i]= 1.0 - (Sim.APAR_Tree[i]/Met_c.PAR[i])
    if abs(Sim.Transmittance_Tree[i])==Inf
        Sim.Transmittance_Tree[i]= 1.0
    end
end

function light_model_coffee!(Sim,Parameters,Met_c,i)
    Sim.K_Dif[i]= Parameters.k_Dif
    Sim.K_Dir[i]= Parameters.k_Dir

    #APAR coffee
    Sim.PAR_Trans_Tree[i]= Met_c.PAR[i] - Sim.APAR_Tree[i] # PAR above coffee layer
    Sim.APAR_Dif[i]= max(0.0, (Sim.PAR_Trans_Tree[i] * Met_c.FDiff[i]) * (1.0 - exp(-Sim.K_Dif[i] * Sim.LAI[i])))
    APAR_Dir= max(0.0,(Sim.PAR_Trans_Tree[i] * (1.0 - Met_c.FDiff[i])) * (1.0 - exp(-Sim.K_Dir[i] * Sim.LAI[i])))
    # APAR_Dir is not part of Sim because it can be easily computed by Met_c.PARm2d1-Sim.APAR_Dif
    Sim.APAR[i]= APAR_Dir + Sim.APAR_Dif[i]
    Sim.PAR_Trans[i]= Sim.PAR_Trans_Tree[i] - Sim.APAR[i] # PAR above soil layer
end

function energy_model_coffee!(Sim,Parameters,Met_c,i)

    # Transpiration Coffee
    Sim.T_Coffee[i]= Base.invokelatest(Parameters.T_Coffee,Sim,Met_c,i)
    # Sensible heat Coffee
    Sim.H_Coffee[i]= Base.invokelatest(Parameters.H_Coffee,Sim,Met_c,i)

    # Tcanopy Coffee : using bulk conductance if no trees, interlayer conductance if trees
      # Source: Van de Griend and Van Boxel 1989.
      if Sim.Height_Tree[i] > Parameters.Height_Coffee
        Sim.Gb_air_canopy[i]= G_interlay(Wind= Met_c.WindSpeed[i], ZHT = Parameters.ZHT, LAI_top= Sim.LAI_Tree[i], LAI_bot= Sim.LAI[i],
                                         Z_top= Sim.Height_Canopy[i], extwind = Parameters.extwind)  
      else
        Sim.G_bulk[i]= G_bulk(Wind = Met_c.WindSpeed[i], ZHT = Parameters.ZHT, Z_top = Sim.Height_Canopy[i],
                              LAI = Sim.LAI[i], extwind = Parameters.extwind)
        Sim.Gb_air_canopy[i]= Sim.G_bulk[i]
      end   

      Sim.Gb_h[i]= Gb_h(Wind = Met_c.WindSpeed[i], wleaf= Parameters.wleaf, LAI_lay=Sim.LAI[i], LAI_abv=Sim.LAI_Tree[i],
                        ZHT = Parameters.ZHT, Z_top = Sim.Height_Canopy[i], extwind= Parameters.extwind)
      
      Sim.TairCanopy[i]= Sim.TairCanopy_Tree[i] + ((Sim.H_Coffee[i] + Sim.H_Soil[i]) * Parameters.MJ_to_W) / 
                         (Sim.air_density_Tree[i] *  Parameters.cp * Sim.Gb_air_canopy[i])
      # NB: if no trees, TairCanopy= Tair (see initialization.jl)

      Sim.air_density[i]= air_density(Sim.TairCanopy[i],Met_c.Pressure[i] / 10.0)

      Sim.Tleaf_Coffee[i]= Sim.TairCanopy[i] + (Sim.H_Coffee[i] * Parameters.MJ_to_W) / 
                            (air_density(Sim.TairCanopy[i],Met_c.Pressure[i] / 10.0) *  Parameters.cp  * Sim.Gb_h[i])

      # Recomputing soil temperature knowing TairCanopy
      Sim.TSoil[i]= Sim.TairCanopy[i]+(Sim.H_Soil[i] * Parameters.MJ_to_W) / 
        (Sim.air_density[i] *  Parameters.cp * 
           G_soilcan(Wind= Met_c.WindSpeed[i], ZHT=Parameters.ZHT, Z_top= max(Sim.Height_Tree[i], Parameters.Height_Coffee),
                     LAI = Sim.LAI_Tree[i] + Sim.LAI[i], extwind= Parameters.extwind))
  
      Sim.DegreeDays_Tcan[i]= GDD(Sim.TairCanopy[i], Parameters.MinTT, Parameters.MaxTT)
end


function energy_model_tree!(Sim,Parameters,Met_c,i)
   
    # Transpiration Coffee
    Sim.T_Tree[i]= Base.invokelatest(Parameters.T_Tree,Sim,Met_c,i)
    # Sensible heat Coffee
    Sim.H_Tree[i]= Base.invokelatest(Parameters.H_Tree,Sim,Met_c,i)

    Sim.G_bulk[i]= G_bulk(Wind= Met_c.WindSpeed[i], ZHT= Parameters.ZHT,
                          LAI= Sim.LAI_Tree[i], extwind= Parameters.extwind,
                          Z_top= Sim.Height_Tree[previous_i(i)])

    # Computing the air temperature in the shade tree layer:
    Sim.TairCanopy_Tree[i]= 
      Met_c.Tair[i] + (Sim.H_Tree[i] * Parameters.MJ_to_W) /
      (air_density(Met_c.Tair[i],Met_c.Pressure[i] / 10.0) * Parameters.cp * Sim.G_bulk[i])
    # NB : using WindSpeed because wind extinction is already computed in G_bulk (until top of canopy).
  
    Sim.Gb_h_Tree[i]= Gb_h(Wind = Met_c.WindSpeed[i], wleaf= Parameters.wleaf_Tree, LAI_lay= Sim.LAI_Tree[i],
                        LAI_abv= 0,ZHT = Parameters.ZHT, Z_top = Sim.Height_Tree[previous_i(i)],
                        extwind= Parameters.extwind)

    Sim.Tleaf_Tree[i]=
      Sim.TairCanopy_Tree[i] + (Sim.H_Tree[i]*Parameters.MJ_to_W) /
      (air_density(Met_c.Tair[i],Met_c.Pressure[i] / 10.0) * Parameters.cp * Sim.Gb_h_Tree[i])        

    Sim.air_density_Tree[i]= air_density(Sim.TairCanopy_Tree[i], Met_c.Pressure[i] / 10.0)
end
