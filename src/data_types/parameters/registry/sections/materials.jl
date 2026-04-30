function get_materials_parameters()::NamedTuple
    return @params (
        # *****************************
        # Materials Parameters
        # *****************************
        
        # Materials - Stress limits powerlaw parameters
        powerlaw_stress_min = ParameterFloat(
            1e4, "Pa", "Minimum stress limit for power law"),

        # Materials - Hydrothermal parameters
        iuse_hydrothermal = ParameterInt(
            0, "None", 
            "Hydrothermal circulation: 0 off; 1 on (iuse_topo must be 1)"
            ),
        iuse_melt_lens = ParameterInt(
            0, "None", "Include shallow melt lens effect: 0 off; 1 on"),
        hydrothermal_smoothing_factor = ParameterFloat(
            0.75, "None", "Hydrothermal smoothing factor"),
        hydrothermal_nusselt_number = ParameterFloat(
            4.0, "None", "Hydrothermal Nusselt number"),
        hydrothermal_max_temperature = ParameterFloat(
            600.0, "C", "Hydrothermal max temperature in Celsius"),
        hydrothermal_max_depth = ParameterFloat(
            6000.0, "m", "Hydrothermal max depth in meters"),
        iuse_plastic_strain_rate_for_hydrothermal = ParameterInt(
            0, "None", 
            "Use plastic strain rate for hydrothermal circulation: 0 off; 1 on"
            ),
        hydrothermal_decay_length = ParameterFloat(
            50000.0, "m", "Hydrothermal decay length in meters"),
        hydrothermal_buffer_distance = ParameterFloat(
            0.0, "m", "Hydrothermal buffer distance in meters"),
        hydrothermal_plastic_strain_rate_reference = ParameterFloat(
            1e-13, "1/s", 
            "Hydrothermal plastic strain rate reference in 1/s"
            ),
        iuse_plastic_strain_for_hydrothermal = ParameterInt(
            0, "None", 
            "Use plastic strain for hydrothermal circulation: 0 off; 1 on"
            ),
        hydrothermal_plastic_strain_reference = ParameterFloat(
            0.5, "None", 
            "Hydrothermal plastic strain reference"
            ),
        sediment_thickness_threshold = ParameterFloat(
            2500.0, "m", "Sediment thickness threshold in meters"),

        # Materials - Melt damage parameters
        iuse_melt_damage = ParameterInt(
            0, "None", 
            "Integer flag that activates melt damage model where the damage "
            *"factor is calculated for markers in the melt damage zone above "
            *"partially molten mantle."
            ),
        melt_damage_distance = ParameterFloat(
            5000.0, "meters", 
            "Distance in meters from the shallow partially molten mantle where the damage "
            *"factor is calculated for markers in the melt damage zone."
            ),
        melt_damage_factor = ParameterFloat(
            10.0, "None", 
            "Maximum melt damage factor. Friction angle is divided by this factor "
            *"to account for melt damage."
            ),
        melt_damage_taper_distance = ParameterFloat(
            1000.0, "meters", 
            "Distance in meters over which melt damage factor tapers to one."
            ),
        iuse_probabilistic_melt_damage = ParameterInt(
            0, "None", 
            "Integer flag that activates probabilistic melt damage model where "
            *"the damage factor is calculated probabilistically for markers in "
            *"the melt damage zone."
            ),
        maximum_damage_probability = ParameterFloat(
            0.75, "fraction", 
            "Maximum probability (fraction) of melt damage for markers in the central melt damage zone."
            ),
        intermediate_damage_probability = ParameterFloat(
            0.3, "fraction", 
            "Intermediate probability (fraction) of melt damage for markers in the central melt damage zone."
            ),
        magmatic_crust_height_threshold = ParameterFloat(
            1000.0, "meters", 
            "Threshold height in meters for magmatic crust used in linear melt damage probability model"
            ),
        magmatic_crust_height_minimum = ParameterFloat(
            3000.0, "meters", 
            "Minimum height in meters for magmatic crust used in linear melt damage probability model"
            ),
        magmatic_crust_height_maximum = ParameterFloat(
            10000.0, "meters", 
            "Maximum height in meters for magmatic crust used in linear melt damage probability model."
            ),
        magmatic_crust_height_intermediate = ParameterFloat(
            2000.0, "meters", 
            "Intermediate height in meters for magmatic crust used in linear melt damage probability model."
            ),
        density_dike_fluid = ParameterFloat(
            2750.0, "kg/m^3",
            "Reference density in kg/m^3 for dike fluid in the melt-damage bulk density correction "
            *"(Birch-Murnaghan EOS at marker pressure)."
            ),
        dike_fluid_marker_fraction = ParameterFloat(
            0.0, "fraction",
            "Volume fraction of dike fluid blended into bulk density when melt damage is active."
            ),

        # Materials - Serpentinization parameters
        iuse_serpentinization = ParameterInt(
            0, "None", 
            "Serpentinization model: 0 off; 1 on. Topography must be activated for serpentinization "
            *"(i.e. `iuse_topo` must be equal to 1). See [EarthBox.SurfaceProcesses.Topography.initialize!](@ref)) "
            *"for more details on topography initialization."
            ),
        serpentinization_temperature = ParameterFloat(
            600.0, "C", 
            "Temperature in Celsius that controls when serpentinization rate reaches maximum and then rapidly "
            *"decreases. Typical values from the literature are around 613.15K (340.15C) "
            *"(e.g. [mezri24](@citet))."
            ),
        maximum_serpentinization_depth = ParameterFloat(
            20000.0, "m", 
            "Maximum submud serpentinization depth in meters"
            ),
        maximum_serpentinization_rate = ParameterFloat(
            1e-10, "1/s", 
            "Maximum serpentinization rate (1/s). Typical values from the literature are around "
            *"1e-10 to 1e-11 1/s (e.g. [mezri24](@citet))."
            ),
        nominal_strain_rate_serpentinization = ParameterFloat(
            1e-13, "1/s", 
            "Nominal plastic strain rate for serpentinization (1/s). Plastic strain rate beyond "
            *"which the effect of plastic strain rate on serpentinization rate rapidly increases "
            *"and then becomes constant (e.g. Merzi et al., 2024). Typical values from the "
            *"literature are around 1e-12 to 1e-13 1/s (e.g. [mezri24](@citet))."
            ),
        # Materials - Description parameters
        itype_mat = ParameterInt(
            0, "None", "Material and marker initialization model option"),
        stype_mat = ParameterStr(
            "None", "None", "Material and marker initialization model option name"),
        itype_plasticity = ParameterInt(
            0, "None", "Plasticity model type: 0 = viscoelastic; 1 = purely elastic"),
        stype_plasticity = ParameterStr(
            "None", "None", "Plasticity model name"),

        # Materials - Random friction parameters
        iuse_random_fric = ParameterInt(
            0, "None", 
            "Randomize initial friction coefficient using delta_fric_coef: 0 off; 1 on"
            ),
        delta_fric_coef = ParameterFloat(
            0.1, "None", 
            "Friction coefficient perturbation used to randomize initial marker friction coefficients"
        ),
        central_delta_fric_coef = ParameterFloat(
            0.1, "None",
            "Friction coefficient perturbation used in the central weakening model "
            *"for randomizing initial marker friction coefficients"
        ),
        central_weakening_probability = ParameterFloat(
            0.5, "None",
            "Maximum weakening probability at the center of the central weakening zone (0 to 1)"
        ),
        iuse_random_fric_time = ParameterInt(
            0, "None", 
            "Randomize marker friction coefficients with each time step: 0 off; 1 on"
            ),
        randomization_factor = ParameterFloat(
            10.0, "None", 
            "Randomization factor used to randomize marker friction coefficients with each time step"
            ),

        # Materials - Stress limits yield parameters
        yield_stress_min = ParameterFloat(
            0.0, "Pa", "Minimum stress limit for plastic failure in Pa"),
        yield_stress_max = ParameterFloat(
            1e38, "Pa", "Maximum stress limit for plastic failure in Pa"),
        iuse_fluid_pressure_for_yield = ParameterInt(
            0, "None", 
            "Flag to use fluid pressure for yield stress calculation"
            ),
        plastic_healing_rate = ParameterFloat(
            0.0, "1/s", "Rate of plastic healing in 1/s"),

        # Materials - Boundary friction parameters
        boundary_friction_width = ParameterFloat(
            0.0, "m", "Width of boundary friction zone"),
        boundary_friction_angle = ParameterFloat(
            0.0, "Degrees", "Friction angle of boundary friction zone"),
        boundary_cohesion = ParameterFloat(
            0.0, "Pa", "Cohesion of boundary friction zone"),

        # Materials - Compaction parameters
        iuse_sed_porosity = ParameterInt(
            0, "None", "Depth-dependent sediment porosity: 0 off; 1 on"),
        conductivity_water = ParameterFloat(
            0.61, "W/m/K", "Thermal conductivity of water"),
        density_water = ParameterFloat(
            1000.0, "kg/m^3", "Density of water"),
        heat_capacity_water = ParameterFloat(
            3000.0, "J/K/kg", "Heat capacity of water"),
        porosity_at_mudline = ParameterFloat(
            0.0, "fraction", "Initial porosity at the sediment-water interface"),

        # Materials - Softening parameters
        iuse_viscous_strain_soft = ParameterInt(
            0, "None", "Viscous strain softening: 0 off; 1 on"),
        vsoftfac = ParameterFloat(
            30.0, "None", "Softening factor multiplied by power-law pre-exponential"),

        # Materials - Viscosity limits parameters
        viscosity_min = ParameterFloat(
            1e18, "Pa.s", "Minimum viscosity"),
        viscosity_max = ParameterFloat(
            1e25, "Pa.s", "Maximum viscosity"),

        )
end
