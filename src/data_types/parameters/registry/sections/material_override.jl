function get_material_override_parameters()::NamedTuple
    return @params (
        #******************
        # Material Override
        # *****************
        # Material override parameters are used to override material properties that are defined in the 
        # material library file.

        # Material Override: Friction Angles
        friction_angle_initial_strong_zone = ParameterFloat(
            0.0, "degrees", 
            "Initial friction angle for strong zones"
            ),
        friction_angle_final_strong_zone = ParameterFloat(
            0.0, "degrees", 
            "Final friction angle for strong zones"
            ),
        friction_angle_initial_solidified_basalt = ParameterFloat(
            0.0, "degrees", 
            "Initial friction angle for solidified basalt"
            ),
        friction_angle_final_solidified_basalt = ParameterFloat(
            0.0, "degrees", 
            "Final friction angle for solidified basalt"
            ),
        friction_angle_initial_oceanic_crust = ParameterFloat(
            0.0, "degrees", 
            "Initial friction angle for oceanic crust"
            ),
        friction_angle_final_oceanic_crust = ParameterFloat(
            0.0, "degrees", 
            "Final friction angle for oceanic crust"
            ),
        friction_angle_initial_sediment = ParameterFloat(
            0.0, "degrees", 
            "Initial friction angle for sediment"
            ),
        friction_angle_final_sediment = ParameterFloat(
            0.0, "degrees", 
            "Final friction angle for sediment"
            ),

        # Material Override: Latent Heat and Solidus and Liquidus
        latent_heat_mantle = ParameterFloat(
            400000.0, "J/kg", 
            "Latent heat for mantle"
            ),
        latent_heat_oceanic_crust = ParameterFloat(
            400000.0, "J/kg", 
            "Latent heat for oceanic crust"
            ),
        mantle_solidus = ParameterStr(
            "PeridotiteKatz2003", "None", 
            "String name of solidus model for mantle"
            ),
        mantle_liquidus = ParameterStr(
            "PeridotiteKatz2003", "None", 
            "String name of liquidus model for mantle"
            ),
        layered_gabbro_solidus = ParameterStr(
            "GabbroGlacier", "None", 
            "String name of solidus model for layered gabbro"
            ),
        layered_gabbro_liquidus = ParameterStr(
            "GabbroGlacier", "None", 
            "String name of liquidus model for layered gabbro"
            ),
        gabbro_solidus = ParameterStr(
            "GabbroGerya2010", "None", 
            "String name of solidus model for gabbro"
            ),
        gabbro_liquidus = ParameterStr(
            "GabbroGerya2010", "None", 
            "String name of liquidus model for gabbro"
            ),

        # Material Override: Pre-exponential factors
        scale_factor_crustal_dislocation_creep = ParameterFloat(
            1.0, "None", 
            "Scale factor for pre-exponential term for crustal dislocation creep"
            ),
        scale_factor_mantle_dislocation_creep = ParameterFloat(
            1.0, "None", 
            "Scale factor for pre-exponential term for mantle dislocation creep"
            ),
        scale_factor_mantle_diffusion_creep = ParameterFloat(
            1.0, "None", 
            "Scale factor for pre-exponential term for mantle diffusion creep"
            ),
        scale_factor_oceanic_crust_dislocation_creep = ParameterFloat(
            1.0, "None", 
            "Scale factor for pre-exponential term for oceanic crust dislocation creep"
            ),
        scale_factor_solidified_basalt_dislocation_creep = ParameterFloat(
            1.0, "None", 
            "Scale factor for pre-exponential term for solidified basalt dislocation creep"
            ),
        scale_factor_sediment_dislocation_creep = ParameterFloat(
            1.0, "None", 
            "Scale factor for pre-exponential term for sediment dislocation creep"
            ),

        # Material Override: Strain limits for weakening
        strain_initial_strong_zone = ParameterFloat(
            0.0, "None", 
            "Initial reference strain used for linear strain weakening model for strong zones"
            ),
        strain_final_strong_zone = ParameterFloat(
            0.0, "None", 
            "Final reference strain used for linear strain weakening model for strong zones"
            ),
        strain_initial_basalt = ParameterFloat(
            0.0, "None", 
            "Initial reference strain used for linear strain weakening model for solidified basalt"
            ),
        strain_final_solidified_basalt = ParameterFloat(
            0.0, "None", 
            "Final reference strain used for linear strain weakening model for solidified basalt"
            ),
        strain_initial_oceanic_crust = ParameterFloat(
            0.0, "None", 
            "Initial reference strain used for linear strain weakening model for oceanic crust"
            ),
        strain_final_oceanic_crust = ParameterFloat(
            0.0, "None", 
            "Final reference strain used for linear strain weakening model for oceanic crust"
            ),
        strain_initial_sediment = ParameterFloat(
            0.0, "None", 
            "Initial reference strain used for linear strain weakening model for sediment"
            ),
        strain_final_sediment = ParameterFloat(
            0.0, "None", 
            "Final reference strain used for linear strain weakening model for sediment"
            ),

        )
end
