function get_topography_parameters()::NamedTuple
    return @params (
        # *****************************
        # Topography Parameters
        # *****************************
        
        # Topography - Reference lithosphere parameters
        thickness_upper_continental_crust_ref = ParameterFloat(
            22000.0, "m", 
            "Thickness of upper continental crust of the reference lithosphere in meters"
            ),
        thickness_lower_continental_crust_ref = ParameterFloat(
            10000.0, "m", 
            "Thickness of lower continental crust of the reference lithosphere in meters"
            ),
        thickness_lithosphere_ref = ParameterFloat(
            125_000.0, "m", 
            "Thickness of lithosphere for the reference lithosphere in meters"
            ),
        gridy_spacing_ref = ParameterFloat(
            100.0, "m", 
            "Spacing of grid cells in y-direction for the reference lithosphere in meters"
            ),
        temperature_top_ref = ParameterFloat(
            0.0, "C", 
            "Temperature at top of lithosphere for the reference lithosphere in Celsius"
            ),
        temperature_moho_ref = ParameterFloat(
            600.0, "C", 
            "Temperature at Moho of the reference lithosphere in Celsius"
            ),
        temperature_base_lith_ref = ParameterFloat(
            1330.0, "C", 
            "Temperature at base of lithosphere for the reference lithosphere in Celsius"
            ),
        adiabatic_gradient_ref = ParameterFloat(
            0.4, "K/km", 
            "Adiabatic gradient below the reference lithosphere in K/km"
            ),
        iuse_linear_segments = ParameterInt(
            1, "None", 
            "Integer flag that controls the type of temperature profile used in "
            *"the reference lithosphere model. If `iuse_linear_segments = 1`, a temperature profile "
            *"with four linear segments is used that is controlled by the following user defined "
            *"parameters: `temperature_top_ref`, `temperature_moho_ref`, and `temperature_base_lith_ref`. "
            *"If `iuse_linear_segments = 0`, the temperature profile is calculated using an analytical 3-layer model "
            *"defined during marker temperature. See the `AnalyticalThreeLayer` model description at "
            *"[MarkerTemperature](@ref)). Default is 1."
            ),

        # Topography - Downhill diffusion parameters
        downhill_diff_elev_max = ParameterFloat(
            1.0, "m", "d2Yt - max elevation"),
        transport_length = ParameterFloat(
            1.0, "m", "Transport length scale in meters"),
        topo_diff_coef = ParameterFloat(
            0.0, "m^2/s", "Topography diffusion coefficient in m^2/s"),
        subaerial_slope_diffusivity = ParameterFloat(
            0.0, "m^2/s", 
            "Subaerial slope diffusivity in m^2/s. Typical value used in the literature is 0.25 "
            *"m^2/yr (7.9e-9 m^2/s) (e.g. Andres-Martinez et al., 2019; Armitage et al., 2015)"
            ),
        precipitation_rate = ParameterFloat(0.0, "m/s", 
            "Precipitation rate in m/s used to calculate water flux in drainage basins for "
            *"fluvial transport diffusivity model. Precipitation rate (m/s). Used to "
            *"calculate water flux in drainage basins. Typical value used in the literature is "
            *"1 m/yr (3.2e-8 m/s) (e.g. Andres-Martinez et al.; 2019; Huffman et al., 2009)."
            ),
        subaerial_transport_coefficient = ParameterFloat(
            0.0, "None", 
            "Subaerial discharge transport coefficient. Used to calculate effective subaerial fluvial "
            *"transport diffusivity that includes slope diffusivity, precipitation rate and downstream "
            *"distances. Typical values used in the literature are 1e-4 (low transport) to 1e-2 (high "
            *"transport) (e.g. Andres-Martinez et al., 2019; Armitage et al., 2015)."
            ),
        submarine_slope_diffusivity = ParameterFloat(
            0.0, "m^2/s", 
            "Maximum submarine slope diffusivity in m^2/s used in diffusivity model that decays exponentially "
            *"with water depth. Typical value used in the literature is 100 m^2/yr (3.2e-9 m^2/s) (e.g. "
            *"Andres-Martinez et al., 2019; Kaufman et al., 1991)."
            ),
        submarine_diffusion_decay_depth = ParameterFloat(
            0.0, "m", 
            "Submarine diffusion decay depth (m). Typical value used in the "
            *"literature is 1000-2000 m (e.g. Andres-Martinez et al., 2019; Kaufman"
            *"et al., 1991; Perez-Gussinye et al., 2020)."
            ),
        number_of_transport_timesteps_per_model_timestep = ParameterInt(
            5, "None", 
            "Number of transport time steps per model time step."
            ),
        transport_timestep = ParameterFloat(
            NaN, "s", 
            "Transport timestep (s). Typical value used in the literature is "
            *"1000 years (3.1536e10 s) (e.g. Andres-Martinez et al., 2019)."
            ),
        iuse_compaction_correction = ParameterInt(
            0, "None", 
            "Integer flag activating compaction correction during sediment transport where newly "
            *"deposited sediment thickness and pre-existing sediment thickness "
            *"are corrected for compaction."
            ),

        # Topography - Depo and erosion rates parameters
        erosion_rate = ParameterFloat(
            0.0, "m/s", 
            "Erosion rate above water level in m/s"
            ),
        sedimentation_rate = ParameterFloat(
            0.0, "m/s", 
            "Sedimentation rate below water level in m/s"
            ),
        pelagic_sedimentation_rate = ParameterFloat(
            0.0, "m/s", 
            "Pelagic sedimentation rate (m/s). Typical value used in the "
            *"literature are 0.3 mm/yr (syn-rift) to 0.01 mm/yr (post-rift) "
            *"(e.g. Perez-Gussinye et al., 2020)"
            ),
        pelagic_sedimentation_rate_reduction_factor = ParameterFloat(
            1.0, "None", 
            "The pelagic sedimentation rate is divided by this factor after "
            *"the specified pelagic sedimentation rate reduction time."
            ),
        pelagic_sedimentation_rate_reduction_time = ParameterFloat(
            5000.0, "Myr", 
            "Time in Myr after which the pelagic sedimentation rate is reduced "
            *"by the pelagic_sedimentation_rate_reduction_factor."
            ),
        salt_deposition_rate = ParameterFloat(
            0.0, "m/s", "Salt deposition rate in m/s"),
        salt_start_time = ParameterFloat(
            1000.0, "Myr", "Salt start time in Myr"),
        salt_end_time = ParameterFloat(
            1001.0, "Myr", "Salt end time in Myr"),

        # Topography - Sea level parameters
        itype_sealevel = ParameterInt(
            0, "None", "Integer flag controlling the sealevel model"),
        stype_sealevel = ParameterStr(
            "None", "None", "Sealevel option name"),
        y_water_ini = ParameterFloat(
            10000.0, "m", "Initial and current seawater level in meters"),
        base_level_shift = ParameterFloat(
            0.0, "m", 
            "Shift in base level in meters. Positive values shift the base level down in the +y direction."
            ),
        base_level_shift_end_time = ParameterFloat(
            -1.0, "Myr", "End time for base level shift in Myr"),
        y_sealevel = ParameterFloat(
            0.0, "m", "Seawater level y-coordinate in meters"),

        # Topography - Model options parameters
        iuse_topo = ParameterInt(
            0, "None", 
            "Integer flag used to enable topography model: 0 = off; 1 = on"
            ),
        itype_topo = ParameterInt(
            0, "None", 
            "ID for topography node advection option: 0 = Runge Kutta with interpolation; 1 = antidiffusion"
            ),
        stype_topo = ParameterStr(
            "None", "None", 
            "Topography advection option name"
            ),
        iuse_downhill_diffusion = ParameterInt(
            0, "None", 
            "Integer flag that activates sediment transport model with downhill diffusion "
            *"(0 = off, 1 = on). If option is set to off uniform erosion and deposition rates will be used"
            ),
        iuse_salt_deposition = ParameterInt(
            0, "None",
            "Integer flag that activates salt deposition model: 0 = off; 1 = on"
            ),

        # Topography - Topo grid parameters
        topo_xsize = ParameterFloat(
            100000.0, "m", "Topography model size in horizontal direction in meters"),
        toponum = ParameterInt(
            1001, "None", "Number of nodes in the x-direction of the topography model"),
        dx_topo = ParameterFloat(
            NaN, "m", "Topography grid step size in meters"),
        nsmooth_top_bottom = ParameterInt(
            2, "None", 
            "Number of topography nodes that will be used to calculate a running average for the "
            *"calculated y-coordinate of the top and bottom of marker layers"
            ),
        marker_search_factor = ParameterFloat(
            2.0, "None", 
            "Search factor used to calculate the search radius for determining the shallowest and "
            *"deepest markers associated with a given topography node"
            ),

        )
end
