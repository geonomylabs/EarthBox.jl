function get_geometry_parameters()::NamedTuple
    return @params (
        # *****************************
        # Geometry parameters
        # *****************************

        # Geometry - Sticky air parameters
        thick_air = ParameterFloat(
            NaN, "m", "Thickness in meters of sticky air in meters"),
        # Geometry - Sandbox parameters
        nsand_layers = ParameterInt(
            0, "None", 
            "Number of sand layers used in sandbox extension and shortening models"
            ),
        y_sand_air_interface = ParameterFloat(
            0.0, "m", 
            "Y-location in meters of interface between sand and air along the left boundary used in "
            *"sandbox extension and shortening models"
            ),
        y_top_microbeads = ParameterFloat(
            0.0, "m", 
            "Y-location in meters of top of microbeads used in sandbox shortening models"
            ),
        y_bottom_microbeads = ParameterFloat(
            0.0, "m", 
            "Y-location in meters of bottom of microbeads used in sandbox shortening models"
            ),
        x_left_ramp = ParameterFloat(
            0.0, "m", 
            "X-location in meters of left-edge of ramp used in sandbox shortening models"
            ),
        x_right_ramp = ParameterFloat(
            0.0, "m", 
            "X-location in meters of right-edge of ramp used in sandbox shortening models"
            ),
        ramp_dip_deg = ParameterFloat(
            0.0, "Degrees", 
            "Dip of ramp in degrees used in sandbox shortening models used in sandbox shortening models"
            ),
        pdms_layer_width = ParameterFloat(
            0.0, "m", 
            "Width in meters of PDMS (polydimethylsiloxane) layer located in the central bottom part "
            *"of model used in sandbox extension models"
            ),
        pdms_layer_thickness = ParameterFloat(
            0.0, "m", 
            "Thickness in meters of PDMS (polydimethylsiloxane) layer used in sandbox extension models"
            ),

        # Geometry - Crustal hole parameters
        xhole_start = ParameterFloat(
            1000.0, "m", "x-coordinate in meters of the start of the crustal hole"),
        xhole_middle = ParameterFloat(
            1100.0, "m", "x-coordinate in meters of the middle of the crustal hole"),
        xhole_end = ParameterFloat(
            1200.0, "m", "x-coordinate in meters of the end of the crustal hole"),
        xhole_depth = ParameterFloat(
            10.0, "m", "Depth in meters of the crustal hole"),

        # Geometry - Mobile wall parameters
        x_left_mobile_wall = ParameterFloat(
            0.0, "m", 
            "X-location in meters of left-edge of mobile wall used in sandbox extension and shortening models"),
        x_right_mobile_wall = ParameterFloat(
            0.0, "m", 
            "X-location in meters of right-edge of mobile wall used in sandbox extension and shortening models"),
        y_top_mobile_wall = ParameterFloat(
            0.0, "m", 
            "Y-location in meters of top-edge of mobile wall used in sandbox extension and shortening models"),
        y_bottom_mobile_wall = ParameterFloat(
            0.0, "m", 
            "Y-location in meters of bottom-edge of mobile wall used in sandbox extension and shortening models"),
        plate_extension_width = ParameterFloat(
            0.0, "m", 
            "Width in meters of plate extension used in sandbox extension models"),
        plate_extension_thickness = ParameterFloat(
            0.0, "m", 
            "Thickness in meters of plate extension used in sandbox extension models"),

        # Geometry - Rayleigh Taylor parameters
        depth_interface_h1 = ParameterFloat(
            100000.0, "m", 
            "Depth in meters of the interface for Rayleigh-Taylor instability"
            ),
        wave_length_lambda = ParameterFloat(
            100000.0, "m", 
            "Wavelength in meters of the perturbation for Rayleigh-Taylor instability"
            ),
        amplitude_initial = ParameterFloat(
            500.0, "m", 
            "Initial amplitude in meters of the perturbation for Rayleigh-Taylor instability"
            ),

        # Geometry - Earth layering parameters
        thick_upper_crust = ParameterFloat(
            NaN, "m", "Thickness in meters of upper crust"),
        thick_lower_crust = ParameterFloat(
            NaN, "m", "Thickness in meters of lower crust"),
        thick_upper_lith = ParameterFloat(
            NaN, "m", "Thickness in meters of upper mantle lithosphere"),
        thick_middle_lith = ParameterFloat(
            NaN, "m", "Thickness in meters of middle lithosphere"),
        thick_lower_lith = ParameterFloat(
            NaN, "m", "Thickness in meters of lower mantle lithosphere"),
        thick_lith = ParameterFloat(
            NaN, "m", "Thickness in meters of lithosphere"),
        thick_crust = ParameterFloat(
            NaN, "m", "Thickness in meters of crust"),
        thick_mantle_lith = ParameterFloat(
            NaN, "m", "Thickness in meters of mantle lithosphere"),

        # Geometry - Fracture zone parameters
        sediment_thickness = ParameterFloat(
            0.0, "m", "Thickness in meters of sediment layer"),
        basaltic_oceanic_crust_thickness = ParameterFloat(
            0.0, "m", 
            "Thickness in meters of basaltic oceanic crust layer"
            ),
        gabbroic_oceanic_crust_thickness = ParameterFloat(
            0.0, "m", 
            "Thickness in meters of gabbroic oceanic crust layer"
            ),
        thickness_of_younger_lithosphere = ParameterFloat(
            0.0, "m", "Thickness in meters of younger lithosphere"),
        thickness_of_older_lithosphere = ParameterFloat(
            0.0, "m", "Thickness in meters of older lithosphere"),
        thickness_of_weak_lithosphere = ParameterFloat(
            0.0, "m", "Thickness in meters of weak lithosphere"),
        x_fracture_zone_start = ParameterFloat(
            0.0, "m", "Start of fracture zone in x-direction in meters"),
        x_fracture_zone_end = ParameterFloat(
            0.0, "m", "End of fracture zone in x-direction in meters"),

        # Geometry - Internal velocity zone parameters
        xindex_vx_internal = ParameterInt(
            -1, "None", "X-index for internal vx velocity"),
        yindex_min_vx_internal = ParameterInt(
            -1, "None", "Minimum y-index for internal vx velocity"),
        yindex_max_vx_internal = ParameterInt(
            -1, "None", "Maximum y-index for internal vx velocity"),
        xindex_vy_internal = ParameterInt(
            -1, "None", "X-index for internal vy velocity"),
        yindex_min_vy_internal = ParameterInt(
            -1, "None", "Minimum y-index for internal vy velocity"),
        yindex_max_vy_internal = ParameterInt(
            -1, "None", "Maximum y-index for internal vy velocity"),
        x_vx_internal = ParameterFloat(
            0.0, "m", "X-coordinate in meters for internal vx velocity"),
        y_min_vx_internal = ParameterFloat(
            0.0, "m", "Minimum y-coordinate in meters for internal vx velocity"),
        y_max_vx_internal = ParameterFloat(
            0.0, "m", "Maximum y-coordinate in meters for internal vx velocity"),
        x_vy_internal = ParameterFloat(
            0.0, "m", "X-coordinate in meters for internal vy velocity"),
        y_min_vy_internal = ParameterFloat(
            0.0, "m", "Minimum y-coordinate in meters for internal vy velocity"),
        y_max_vy_internal = ParameterFloat(
            0.0, "m", "Maximum y-coordinate in meters for internal vy velocity"),

        # Geometry - Litho strong zones parameters
        x_left_strong = ParameterFloat(
            100000.0, "m", 
            "x-coordinate in meters of the left edge of the lithospheric strong zone"
            ),
        x_right_strong = ParameterFloat(
            200000.0, "m", 
            "x-coordinate in meters of the right edge of the lithospheric strong zone"
            ),
        iuse_strong_zones = ParameterInt(
            0, "None", "Use strong zones (0=off, 1=on)"),

        # Geometry - Plume parameters
        plume_radius = ParameterFloat(
            200000.0, "m", "Radius in meters of the mantle plume"),
        plume_center_x = ParameterFloat(
            100000.0, "m", "x-coordinate in meters of the plume center"),
        plume_center_y = ParameterFloat(
            600000.0, "m", "y-coordinate in meters of the plume center"),
        plume_head_thick = ParameterFloat(
            100000.0, "m", "Thickness in meters of the plume head"),
        delta_temperature_plume = ParameterFloat(
            100.0, "K", 
            "Temperature difference of the plume relative to background mantle"
            ),
        iuse_plume = ParameterInt(
            0, "None", "Flag for using plume (0=off, 1=on)"),
        plume_start_time = ParameterFloat(
            0.0, "Myr", "Time in Myr when plume starts"),
        iplume_state = ParameterInt(
            0, "None", "State of the plume: 0 - not initiated; 1 = initiated"),

        # Geometry - Weak fault parameters
        fault_dip_degrees = ParameterFloat(
            0.0, "Degrees", 
            "Dip in degrees of the weak zone used to approximate a fault"
            ),
        fault_thickness = ParameterFloat(
            0.0, "m", 
            "Thickness in meters of the weak zone used to approximate a fault"
            ),
        x_initial_fault = ParameterFloat(
            0.0, "m", 
            "Initial x-coordinate in meters of the weak zone used to approximate a fault"
            ),
        fault_height = ParameterFloat(
            0.0, "m", "Height in meters of the weak zone used to approximate a fault"),
        iuse_weak_fault = ParameterInt(
            0, "None", "Use weak fault (0=off, 1=on)"),

        # Geometry - Weak seed parameters
        w_seed = ParameterFloat(
            1000.0, "m", "Width in meters of weak seed"),
        x_seed = ParameterFloat(
            100000.0, "m", "x-location in meters of weak seed"),
        y_seed = ParameterFloat(
            35000.0, "m", "y-location in meters of weak seed"),
        iuse_weak_seed = ParameterInt(
            0, "None", "Use weak seed (0=off, 1=on)"),
        )
end
