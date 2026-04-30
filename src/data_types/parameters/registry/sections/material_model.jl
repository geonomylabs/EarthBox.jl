function get_material_model_parameters()::NamedTuple
    return @params (
        #********************
        # Material Parameters
        # *******************
        # Material parameters are defined in material input and library files.

        # Material Model Parameters
        name = ParameterStr("dummy", "None", "Name of material"),
        ncolors = ParameterInt(0, "None", "Number of colors"),
        nmats = ParameterInt(
            100, "None", "Maximum number of materials allowed in model"),
        material_model_id = ParameterInt(
            0, "None", "Material model ID"),
        matid = ParameterInt(0, "None", "Material ID"),
        mat_name = ParameterStr("dummy", "None", "Material name"),
        mat_type = ParameterStr("dummy", "None", "Material type"),
        mat_domain = ParameterStr("dummy", "None", "Material domain"),
        description = ParameterStr("dummy", "None", "Material description"),

        # Density Parameters
        standard_density = ParameterFloat(
            0.0, "kg/m^3", "Standard density in kg/m^3"),
        thermal_expansion = ParameterFloat(
            0.0, "1/K", "Thermal expansion coefficient in 1/K"),
        compressibility = ParameterFloat(
            0.0, "1/Pa", "Compressibility in 1/Pa"),
        melt_density = ParameterFloat(0.0, "kg/m^3", "Melt density in kg/m^3"),
        viscosity_iso = ParameterFloat(0.0, "Pa.s", 
            "Isoviscous viscosity in Pa.s"),

        # Flow Law Parameters
        flow_type = ParameterInt(0, "None", "Flow law type"),
        flow_stype = ParameterStr("None", "None", "Flow law type name"),

        # Dislocation Creep Parameters
        pre_exponential_dc = ParameterFloat(
            0.0, "1/s/MPa^n", "Pre-exponential factor for dislocation creep in 1/s/MPa^n"),
        stress_exponent_n_dc = ParameterFloat(
            0.0, "None", "Stress exponent for dislocation creep"),
        activation_energy_dc = ParameterFloat(
            0.0, "kJ/mol", "Activation energy for dislocation creep in kJ/mol"),
        activation_volume_dc = ParameterFloat(
            0.0, "J/MPa/mol", "Activation volume for dislocation creep in J/MPa/mol"),

        # Peierls Creep Parameters
        pre_exponential_pei = ParameterFloat(
            0.0, "s^-m1*MPa^-m2", "Pre-exponential factor for Peierls creep in s^-m1*MPa^-m2"),
        stress_exponent_m1_pei = ParameterFloat(
            0.0, "None", "Stress exponent m1 for Peierls creep"),
        stress_exponent_m2_pei = ParameterFloat(
            0.0, "None", "Stress exponent m2 for Peierls creep"),
        peierls_stress = ParameterFloat(
            0.0, "MPa", "Peierls stress in MPa"),

        # Diffusion Creep Parameters
        pre_exponential_difc = ParameterFloat(
            0.0, "1/s/MPa^n", "Pre-exponential factor for diffusion creep in 1/s/MPa^n"),
        activation_energy_difc = ParameterFloat(
            0.0, "kJ/mol", "Activation energy for diffusion creep in kJ/mol"),
        activation_volume_difc = ParameterFloat(
            0.0, "J/MPa/mol", "Activation volume for diffusion creep in J/MPa/mol"),

        # Temperature Dependent Viscosity Parameters
        pre_exponential_td = ParameterFloat(
            0.0, "Pa.s", 
            "Pre-exponential factor for temperature dependent viscosity in Pa.s"
            ),
        activation_energy_td = ParameterFloat(
            0.0, "kJ/mol", 
            "Activation energy for temperature dependent viscosity in kJ/mol"
            ),

        # Blankenbach89 Viscosity Parameters
        viscosity_ref_blankenbach89 = ParameterFloat(
            0.0, "Pa.s", "Reference viscosity for Blankenbach89 in Pa.s"),
        b_term_blankenbach89 = ParameterFloat(
            0.0, "None", "b term for Blankenbach89"),
        c_term_blankenbach89 = ParameterFloat(
            0.0, "None", "c term for Blankenbach89"),

        # Shear Modulus Parameters
        shear_modulus = ParameterFloat(
            0.0, "Pa", "Shear modulus in Pa"),

        # Plasticity Parameters
        cohesion_initial = ParameterFloat(
            0.0, "Pa", "Initial cohesion in Pa"),
        cohesion_final = ParameterFloat(
            0.0, "Pa", "Final cohesion in Pa"),
        friction_angle_initial = ParameterFloat(
            0.0, "Degrees", "Initial friction angle in degrees"),
        friction_angle_final = ParameterFloat(
            0.0, "Degrees", "Final friction angle in degrees"),
        sine_friction_angle_initial = ParameterFloat(
            0.0, "None", "Initial sine friction angle"),
        sine_friction_angle_final = ParameterFloat(
            0.0, "None", "Final sine friction angle"),
        strain_initial = ParameterFloat(
            0.0, "None", "Initial strain"),
        strain_final = ParameterFloat(
            0.0, "None", "Final strain"),

        # Heat Capacity Parameters
        heat_capacity = ParameterFloat(
            0.0, "J/kg/K", "Heat capacity in J/kg/K"),

        # Thermal Conductivity Parameters
        thermal_conductivity_ref = ParameterFloat(
            0.0, "W/m/K", "Reference thermal conductivity in W/m/K"),
        thermal_conductivity_a = ParameterFloat(
            0.0, "W/m", "Thermal conductivity parameter a in W/m"),

        # Radiogenic Heat Production Parameters
        radiogenic_heat_production = ParameterFloat(
            0.0, "W/m^3", "Radiogenic heat production in W/m^3"),

        # Material Color Parameters
        red_fraction = ParameterFloat(
            0.0, "None", "Red color fraction"),
        green_fraction = ParameterFloat(
            0.0, "None", "Green color fraction"),
        blue_fraction = ParameterFloat(
            0.0, "None", "Blue color fraction"),
        
        # Dilatation Angle Parameters
        dilatation_angle = ParameterFloat(
            0.0, "Degrees", "Dilatation angle in degrees"),

        # Melting Parameters
        itype_solidus = ParameterInt(
            -1, "None", "Solidus type ID"),
        stype_solidus = ParameterStr(
            "None", "None", "Solidus type name"),
        itype_liquidus = ParameterInt(
            -1, "None", "Liquidus type ID"),
        stype_liquidus = ParameterStr(
            "None", "None", "Liquidus type name"),
        latent_heat = ParameterFloat(
            400000.0, "J/kg", "Latent heat in J/kg"),

        # Compaction Parameters
        porosity_initial = ParameterFloat(
            0.0, "fraction", "Initial porosity in fraction"),
        porosity_decay_depth = ParameterFloat(
            2001.0, "m", "Porosity decay depth in m"),

        )
end
