function get_boundary_condition_parameters()::NamedTuple
    return @params (
        # BCs - Pressure
        pressure_bc = ParameterFloat(
            10000.0, "Pa", "Pressure boundary condition in Pa."),

        # BCs - Velocity parameters
        velocity = ParameterFloat(
            0.0, "m/s", "General velocity boundary condition"),
        full_velocity_extension = ParameterFloat(
            0.0, "m/s", "Full extension rate in m/s"),
        full_velocity_extension_step1 = ParameterFloat(
            NaN, "m/s", "Full extension rate at step 1 in m/s"),
        full_velocity_extension_step2 = ParameterFloat(
            NaN, "m/s", "Full extension rate at step 2 in m/s"),
        full_velocity_contraction = ParameterFloat(
            0.0, "m/s", "Full contraction rate in m/s"),
        velocity_shear = ParameterFloat(
            0.0, "m/s", "Shear velocity in m/s"),
        velocity_rotation = ParameterFloat(
            0.0, "m/s", "Rotation velocity in m/s"),
        iuse_strain_rate = ParameterInt(
            0, "None", 
            "0 use full extension velocity to define velocity bc's; 1 use strain rate"
            ),
        strain_rate_bc = ParameterFloat(
            1e-15, "1/s", "Strain rate boundary condition in 1/s"),
        vyu = ParameterFloat(
            0.0, "m/s", "Velocity in y direction at upper boundary in m/s"),
        vyl = ParameterFloat(
            0.0, "m/s", "Velocity in y direction at lower boundary in m/s"),
        velocity_internal_x = ParameterFloat(
            0.0, "m/s", "X component of internal velocity in m/s"),
        velocity_internal_y = ParameterFloat(
            0.0, "m/s", "Y component of internal velocity in m/s"),
        plate_thickness = ParameterFloat(
            100000.0, "m",
            "Thickness of plate in meters used for inflow-outflow depth-dependent extension and "
            *"contraction boundary conditions"
            ),
        smoothing_thickness = ParameterFloat(
            10000.0, "m", 
            "Thickness in meters for velocity smoothing between inflow and outflow boundaries"
            ),
        velocity_inflow_left = ParameterFloat(
            0.0, "m/s", "Inflow velocity at left boundary in m/s"),
        velocity_inflow_right = ParameterFloat(
            0.0, "m/s", "Inflow velocity at right boundary in m/s"),
        velocity_inflow_smooth_avg_left = ParameterFloat(
            0.0, "m/s", 
            "Smoothed average inflow velocity at left boundary in m/s"
            ),
        velocity_inflow_smooth_avg_right = ParameterFloat(
            0.0, "m/s", 
            "Smoothed average inflow velocity at right boundary in m/s"
            ),
        iuse_velocity_stop = ParameterInt(
            0, "None", "Flag for velocity stopping: 0 = off; 1 = on"),
        velocity_stop_time = ParameterFloat(
            0.0, "Myr", "Time at which velocity stops in Myr"),
        ivelocity_stop_counter = ParameterInt(
            0, "None", "Counter for velocity stop: 0 = off; 1 = on"),

        # BCs - Temperature parameters
        temperature_top = ParameterFloat(
            0.0, "K", "Temperature along the top boundary of the model domain in Kelvin"),
        temperature_bottom = ParameterFloat(
            0.0, "K", "Temperature at the base of the model domain in Kelvin"),
        temperature_left = ParameterFloat(
            0.0, "K", "Temperature along the left boundary of the model domain in Kelvin"),
        temperature_right = ParameterFloat(
            0.0, "K", "Temperature along the right boundary of the model domain in Kelvin"),
        iuse_bottom_transient = ParameterInt(
            0, "None", 
            "Flag to use transient temperature at the bottom boundary"
            ),
        temperature_bottom_transient = ParameterFloat(
            NaN, "K", 
            "transient temperature at the bottom boundary"
            ),
        temperature_bottom_original = ParameterFloat(
            -99999.0, "K", 
            "Original temperature at the bottom boundary"
            ),
        start_time_bottom_transient = ParameterFloat(
            0.0, "Myr", 
            "Start time (Myr) for transient temperature at the bottom boundary"
            ),
        end_time_bottom_transient = ParameterFloat(
            0.0, "Myr", 
            "End time (Myr) for transient temperature at the bottom boundary"
            ),
        delta_temperature_transient = ParameterFloat(
            NaN, "deltaK", 
            "Temperature perturbation for initial transient temperature at the bottom boundary"
            ),
        temperature_base_lith_warmer_initial = ParameterFloat(
            NaN, "K", 
            "Temperature at the base of the lithosphere"
            ),
        temperature_bottom_warmer_initial = ParameterFloat(
            NaN, "K", 
            "Initial elevated temperature at the bottom of the model"
            ),
        temperature_bottom_cooler_final = ParameterFloat(
            NaN, "K", 
            "Cooler transient temperature at the bottom of the model"
            ),

        # BCs - Velocity step parameters
        iuse_velocity_step = ParameterInt(
            0, "None", 
            "Increase/decrease full extension velocity by factor at specified time"
            ),
        velocity_step_factor = ParameterFloat(
            NaN, "None", 
            "Factor used to increase/decrease velocity at specified time"
            ),
        timestep_adjustment_factor = ParameterFloat(
            1.0, "None", 
            "Factor used to adjust time step when velocity is increased/decreased"
            ),
        velocity_step_time = ParameterFloat(
            NaN, "Myr", "Time of velocity increase/decrease in Myr"),
        velocity_second_step_time = ParameterFloat(
            NaN, "Myr", "Time of second velocity increase/decrease in Myr"),
        ivelocity_step_counter = ParameterInt(
            0, "None", "Integer counter for velocity step"),
        velocity_second_step_factor = ParameterFloat(
            NaN, "None", 
            "Factor used to increase/decrease velocity at specified time"
            ),
        timestep_second_adjustment_factor = ParameterFloat(
            1.0, "None", 
            "Factor used to adjust time step when velocity is increased/decreased"
            ),

        # BCs - Bc option parameters
        itype_bc = ParameterInt(
            0, "None", "Boundary condition option for Stokes and heat equations"),
        stype_bc = ParameterStr(
            "None", "None", "Boundary condition option name"),
        pressure_bc_mode = ParameterFloat(
            0.0, "None", "Pressure bc: 0 upper left cell 1 top and bottom"),
    )
end
