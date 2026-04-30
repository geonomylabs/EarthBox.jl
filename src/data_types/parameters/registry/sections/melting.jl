function get_melting_parameters()::NamedTuple
    return @params (
        # *****************************
        # Melting Parameters
        # *****************************
        
        # Melting - Rheology parameters
        viscosity_melt = ParameterFloat(
            1e17, "Pa.s", "Viscosity of molten rock in Pa.s"),

        # Melting - Extrusion parameters
        iuse_extrusion = ParameterInt(
            0, "None", 
            "Integer flag that activates melt extrusion model: 0 off; 1 on"
            ),
        extrusion_volume_factor = ParameterFloat(
            0.2, "None",
            "Minimum value for the extrusion volume factor. "  
            *"The volume of lava produced at the surface of the model is calculated "
            *"as the product of the extrusion volume factor and the total volume "
            *"of magma produced during a time step. "
            *"The extrusion volume factor changes as a function of the characteristic "
            *"magmatic crust height. The characteristic magmatic crust height is "
            *"the height of new crust that would be formed if all "
            *"the melt in the mantle was extracted and emplaced at the top of the "
            *"model in a column with a width equal to the full extension velocity "
            *"times the time step. "
            *"This minimum value is used when the magmatic crust height is below the minimum "
            *"characteristic magmatic crust height. A linear model is used to define "
            *"the extrusion volume factor when the calculated magmatic crust height is between "
            *"the minimum and maximum characteristic magmatic crust height."
            ),
        extrusion_volume_factor_max = ParameterFloat(
            0.2, "None", 
            "Maximum value for the extrusion volume factor. This maximum "
            *"value is reached when the characteristic magmatic crust height is equal to the "
            *"maximum characteristic magmatic crust height."
            ),
        characteristic_magmatic_crust_height = ParameterFloat(
            0.0, "m", 
            "Characteristic height of magmatic crust in meters calculated based on "
            *"the total volume of magma produced during a given time step and the "
            *"full extensional velocity. The characteristic magmatic crust height "
            *"is the height of new crust that would be formed if all "
            *"the melt in the mantle was extracted and emplaced at the top of the "
            *"model in a column with a width equal to the full extension velocity "
            *"times the time step."
            ),
        characteristic_magmatic_crust_height_min = ParameterFloat(
            8000.0, "m", 
            "Minimum reference value for the characteristic magmatic crust height in meters. When "
            *"the calculated characteristic magmatic crust height is below this value the "
            *"extrusion volume factor is set to the minimum value "
            *"(i.e. `extrusion_volume_factor`)."
            ),
        characteristic_magmatic_crust_height_max = ParameterFloat(
            8000.0, "m", 
            "Maximum reference value for the characteristic magmatic crust height in meters. When "
            *"the calculated characteristic magmatic crust height is above this value the "
            *"extrusion volume factor is set to the maximum value "
            *"(i.e. `extrusion_volume_factor_max`)."
            ),
        characteristic_flow_length_subaerial = ParameterFloat(
            100000.0, "m", 
            "Characteristic flow length for subaerial eruptions in meters. Flow "
            *"distance of basaltic flows on land may be hundreds of kilometers (e.g. Ginkgo "
            *"flow of the Frenchman Springs Member of the Columbia River Basalts "
            *"has a flow length of 500 km)."
            ),
        characteristic_flow_length_submarine = ParameterFloat(
            10000.0, "m", 
            "Characteristic flow length for submarine eruptions in meters."
            ),
        residual_lava_thickness_subaerial = ParameterFloat(
            30.0, "m", 
            "Residual lava thickness for subaerial eruptions in meters"
            ),
        residual_lava_thickness_submarine = ParameterFloat(
            60.0, "m", 
            "Residual lava thickness for submarine eruptions in meters"
            ),
        iuse_random_eruption_location = ParameterInt(
            0, "None", 
            "Integer flag that activates random eruption location: 0 off; 1 on"
            ),
        iuse_normal_eruption_location = ParameterInt(
            0, "None", 
            "Integer flag that activates normal distribution for eruption location: 0 off; 1 on"
            ),
        decimation_factor = ParameterInt(
            1, "None", 
            "Topography grid decimation factor for 1D lava flow model grid"
            ),
        initial_magma_flush_steps = ParameterInt(
            0, "None", 
            "Number of time steps for initial magma flush"
            ),
        magma_flush_factor = ParameterFloat(
            0.8, "None", 
            "Factor that determines the amount of magma to flush for flush during magma flush time steps"
            ),
        width_eruption_domain_fixed = ParameterFloat(
            0.0, "m", 
            "Fixed width of the eruption domain in meters. If this is set to zero, the "
            *"width of the eruption domain is determined by the width of the "
            *"gabbroic molten zone including normal and layered gabbro."
            ),
        width_eruption_domain_fixed_max = ParameterFloat(
            0.0, "m", 
            "Maximum width of the eruption domain. If this is set to zero, the "
            *"width of the eruption domain is determined by the width of the "
            *"gabbroic molten zone including normal and layered gabbro."
            ),
        porosity_initial_lava_flow = ParameterFloat(
            0.0, "fraction", 
            "Initial porosity of lava flow in fraction used to decompact extruded gabbroic "
            *"magma to account for vesicles. This parameter is not used to "
            *"compact lava flows during burial. Instead, compaction properties "
            *"associated with the material type are used. This allows for the "
            *"decompaction of lava flows during eruption to account for porosity "
            *"from vesicles and the preservation of thickness during burial due to "
            *"chemical compaction and rigidity of the basalt matrix."
            ),
        decay_depth_lava_flow = ParameterFloat(
            2000.0, "m", 
            "Depth at which the porosity of lava flow decays to zero in meters. This "
            *"parameter is used to compact lava flows during burial. Instead, "
            *"compaction properties associated with the material type are used. "
            *"This allows for the decompaction of lava flows during eruption to "
            *"account for porosity from vesicles and the preservation of thickness "
            *"during burial due to chemical compaction and rigidity of the basalt "
            *"matrix."
            ),
        time_of_next_eruption_myr = ParameterFloat(
            0.0, "Myr", "Time of next eruption in Myr"),
        eruption_interval_yr = ParameterFloat(
            100000.0, "yr", "Eruption interval in years"),
        iuse_eruption_interval = ParameterInt(
            0, "None", 
            "Integer flag that activates eruption interval: 0 off; 1 on"
            ),

        # Melting - Extraction parameters
        melt_residual = ParameterFloat(
            0.0, "None", "Sum of residual fractional melt particles"),
        ext_vol = ParameterFloat(
            0.0, "m^3", "Volume of extruded material in m^3"),
        xmid_mol = ParameterFloat(
            0.0, "m", "x-location of the molten mantle zone in meters"),
        ytop_mol = ParameterFloat(
            0.0, "m", "y-location of the top of the molten mantle zone in meters"),
        width_mol = ParameterFloat(
            0.0, "m", "Width of the molten mantle zone in meters"),
        ndrainage_basin = ParameterInt(
            1, "None", "Number of drainage basins in the model"),
        smoothing_radius_drainage = ParameterFloat(
            2000.0, "m", 
            "Smoothing radius (meters) used to smooth the top of the partially "
            *"molten mantle domain before calculating local maxima used to define "
            *"migration domains for melt extraction."
            ),
        characteristic_injection_width = ParameterFloat(
            2500.0, "m", 
            "Characteristic width (meters) of magma injection into the shallow "
            *"mantle above the local maximum of the mantle partial melt domain. "
            *"This parameters is used only if use_shallow_mantle_injection is True. "
            *"Note that the actual injection width may be increased if the height "
            *"of injected magma exceeds the magma height limit."
            ),
        magma_height_limit = ParameterFloat(
            30000.0, "m", 
            "Maximum height (meters) of magma that can be injected into the "
            *"mantle domain. This is used to limit the height of magma bodies "
            *"that are injected at the top of the mantle domain by increasing the "
            *"injection width."
            ),
        fractionation_threshold_limit = ParameterFloat(
            2000.0, "m", 
            "Threshold limit (meters) for gabbroic fractionation. This is the "
            *"distance from Moho where the composition of gabbroic magma is changed "
            *"to fractionated or layered gabbroic magma to simulate the effects of "
            *"rapid fractionation in sills and the downward flow of the gabbro "
            *"glacier. This change in composition involves an increase in the solidus "
            *"leading to regions of partially molten fractioned gabbro as opposed "
            *"to large regions of pure gabbroic magma inconsistent with geophysical "
            *"observations."
            ),
        emplacement_temperature = ParameterFloat(
            1473.0, "K", 
            "Temperature at which magma is emplaced in Kelvin. "
            *"This parameter is currently not used in the code. Instead a default value of 1473.0 K is used."
            *"See function `transform_marker_to_magma` for details."
            ),
        number_of_injection_subdomains = ParameterInt(
            10, "None", 
            "The number used to divide the characteristic injection width into "
            *"separate domains that are selected using a normal distribution if "
            *"`iuse_random_injection_subdomain = 1`."
            ),
        maximum_shallow_injection_depth = ParameterFloat(
            200000.0, "m", 
            "Maximum submud depth (meters) for normal shallow mantle injection. "
            *"If the local maximum of the partially molten zone is deeper than "
            *"this value then the injection width is increased by a factor of 5. "
            *"THIS IS CURRENTLY TURNED OFF IN CODE. See function "
            *"extract_partial_melt_and_make_magma_body for details."
            ),
        extraction_fraction = ParameterFloat(
            1.0, "fraction", 
            "Fraction of extracted melt that is extracted from the partially "
            *"molten zone. The remaining melt is left in the partially molten "
            *"zone and refertalizes the mantle."
            ),
        smoothing_radius_fractionation = ParameterFloat(
            10000.0, "m", 
            "Smoothing radius (meters) used to smooth the oceanic Moho for the  "
            *"gabbroic fractionation model whereby fractionation occurs based on "
            *"distance to Moho. Smoothing the Moho avoids spikes associated with "
            *"residual gabbroic particles in the mantle associated with magmatic "
            *"crust thinning."
            ),
        mantle_search_width = ParameterFloat(
            200000.0, "m", 
            "Width (meters) of the mantle domain used to search for the shallowest "
            *"mantle marker in a subset of the model domain to improve computational "
            *"efficiency. This parameter is used only if use_shallow_mantle_injection "
            *"is True and may need to be adjusted depending on the characteristic "
            *"injection width. Adjustments should be made to ensure that the search "
            *"domain is equal to or larger than the injection domain."
            ),
        ndrainage_basin_old = ParameterInt(
            0, "None", "Old number of drainage basins"),
        iuse_melt_compaction = ParameterInt(
            0, "None", 
            "Integer flag that activates melt compaction model in the partially molten zone: 0 off; 1 on."
            ),

        # Melting - Options parameters
        iuse_melting = ParameterInt(
            0, "None", 
            "Integer flag that activates melt model with melt fraction calculation: 0 off; 1 on"
            ),
        iuse_melt_viscosity = ParameterInt(
            0, "None", 
            "Integer flag that activates melt viscosity model where the viscosity of partially "
            *"molten rocks is set to melt viscosity if melt fraction exceeds 10%: 0 off; 1 on"
            ),
        iuse_melt_thermal_props = ParameterInt(
            0, "None", 
            "Integer flag that activates melt to impact thermal properties including density, heat "
            *"capacity in the term rho*cp and thermal expansivity in the adiabatic heating term alpha*T: 0 off; 1 on"
            ),
        iuse_extraction = ParameterInt(
            0, "None", 
            "Integer flag that activates melt extraction model (0 = off, 1 = on). Extraction is required "
            *"for intrusive and extrusive processes."
            ),
        iuse_gabbroic_fractionation = ParameterInt(
            0, "None", 
            "Integer flag that activates gabbroic fractionation model: 0 off; 1 on. This "
            *"model is used to approximate the effects of rapid fractionation in "
            *"sills that produce layered gabbroic magma. This is accomplished by "
            *"transforming gabbroic magma to a fractionated gabbroic magma "
            *"composition once magma enters the crustal domain at the Moho."
            ),
        iuse_shallow_mantle_injection = ParameterInt(
            0, "None", 
            "Integer flag that activates shallow mantle injection model: 0 off; 1 on. This "
            *"model is used to inject mantle-derived magma at the shallowest part "
            *"of the mantle domain. If this option is False then mantle-derived "
            *"magma is injected at the shallowest location of the partially "
            *"molten mantle domain. This option approximates rapid ascent of melt "
            *"from the mantle to the crustal domain via dykes and channelized "
            *"melt networks."
            ),
        iuse_random_injection_subdomain = ParameterInt(
            0, "None", 
            "Integer flag that activates the random injection model (0 = off, 1 = on) whereby "
            *"the characteristic injection width (characteristic_injection_width) is divided into "
            *"subdomains (number_of_injection_subdomains) that are selected using a uniform distribution. "
            *"The shallowest mantle marker within a selected subdomain is converted into magma. If this option "
            *"is False then the shallowest mantle marker within the injection width is converted into magma."
            ),
        iuse_normal_injection_subdomain = ParameterInt(
            0, "None", 
            "Integer flag that activates the normal injection model (0 = off, 1 = on) whereby the "
            *"characteristic injection width (characteristic_injection_width) is divided into "
            *"subdomains (number_of_injection_subdomains) that are selected using a normal distribution. "
            *"The shallowest mantle marker within a selected subdomain is converted into magma. If this option "
            *"is False then a uniform distribution is used to select the subdomain."
            ),
        iuse_depletion_density = ParameterInt(
            0, "None", 
            "Integer flag that activate depletion density model (0 = off, 1 = on) whereby the "
            *"density of rocks that have undergone melt extraction is modified to account for "
            *"the partitioning of heavy elements into the melt."
            ), 
        iuse_exponential_viscosity_reduction = ParameterInt(
            0, "None", 
            "Integer flag that activates exponential viscosity reduction for partially "
            *"molten rocks whereby the effective viscosity of partially molten rocks is "
            *"reduced exponentially with extractable melt fraction for mantle rocks that are impacted by"
            *"melt extraction and melt fraction for crustal rocks that are currently not impacted "
            *"by melt extraction (0 = off, 1 = on). If this parameter is set to 0, the effective "
            *"viscosity of partially molten rocks is set to melt viscosity if melt fraction exceeds 10%."
            ),

        # Melting - Material IDs parameters
        ipmf1 = ParameterInt(
            5, "None", "Material ID of the mantle"),
        ipmf2 = ParameterInt(
            6, "None", "Material ID of the mantle"),
        ipmf3 = ParameterInt(
            7, "None", "Material ID of the mantle"),
        ipmf4 = ParameterInt(
            8, "None", "Material ID of the mantle"),
        ipmf_molten = ParameterInt(
            11, "None", "Material ID of molten mantle"),
        ipmf_solid = ParameterInt(
            12, "None", "Material ID of solidified mantle"),
        ipmf_molten_vol = ParameterInt(
            18, "None", "Material ID of molten extruded mantle"),
        ipmf_solid_vol = ParameterInt(
            19, "None", "Material ID of solidified extruded molten mantle"),

        )
end
