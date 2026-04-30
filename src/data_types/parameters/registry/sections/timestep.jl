function get_timestep_parameters()::NamedTuple
    return @params (
        # *****************************
        # Timestep Parameters
        # *****************************
        
        # Timestep - Boundary displacement stopping parameters
        iuse_boundary_displacement = ParameterInt(0, "None", 
            "Flag to use boundary displacement stopping: 0 = off; 1 = on. When active, the model will "
            *"stop when the boundary displacement reaches the limit specified by `displ_limit`"),
        displ_limit = ParameterFloat(
            500.0, "m", "Displacement limit for stopping the model in meters"),
        iuse_extensional_strain = ParameterInt(
            0, "None",
            "Flag to use extensional strain stopping: 0 = off; 1 = on. When active, the model will "
            *"stop when the extensional strain reaches the limit specified by `strain_limit`"
            ),
        strain_limit = ParameterFloat(
            1.0, "None", "Strain limit for stopping the model"),

        # Timestep - Single increase parameters
        iuse_single_timestep_increase = ParameterInt(
            0, "None", 
            "Flag to use single timestep increase: 0 = off; 1 = on. When active, the time step will be "
            *"increased once during the simulation at the time step specified by ntime_increase"
            ),
        ntime_increase = ParameterInt(100, "None", 
            "Time step number at which to increase the timestep. This is only used if "
            *"`iuse_single_timestep_increase` is 1"
            ),
        ntime_output = ParameterInt(
            100, "None", "Time step number at which to output the model"),
        model_duration_myr = ParameterFloat(
            0.0, "Myr", 
            "Model duration in Myr used to estimate the maximum number of time steps for extensional "
            * "model with symmetric extension. "
            * "If adaptive time stepping is used during the simulation due to displacements that exceed "
            * "the user specified cell fraction then the actual model duration may be shorter than "
            * "`model_duration_myr`. "
            * "A buffer time of 3.0 Myr is added to the model duration when calculating the maximum "
            * "number of time steps assuming that there will be some time step reduction. "
            * "However, the user may need to increase the `model_duration_myr` to achieve the "
            * "desired model duration. If you want full control over the number of time steps then "
            * "set `ntimestep_max` and `timestep_viscoelastic` directly when using `run_time_steps` "
            * "function."
            ),
        ntimestep_max = ParameterInt(
            -9999, "None", "Maximum number of model time steps"),
        timesum = ParameterFloat(
            0.0, "s", "Total model time in seconds in seconds"),
        ntimestep = ParameterInt(
            0, "None", "Model time step counter"),
        timestep_viscoelastic = ParameterFloat(
            NaN, "s", 
            "The viscoelastic time step in seconds used in the viscoelastic formulation of the "
            *"Stokes-continuity equations. The main model time step (`timestep`) is set equal to the "
            *"viscoelastic time step at the start of each time loop iteration prior to defining the "
            *"right-hand side terms that include the time-dependent viscoelastic factor. The "
            *"viscoelastic time step may be modified due to changes in boundary conditions depending "
            *"on user defined options"
            ),
        timestep_viscoelastic_step1 = ParameterFloat(
            NaN, "s", 
            "The viscoelastic time step in seconds for the first velocity step in seconds"
            ),
        timestep_viscoelastic_step2 = ParameterFloat(
            NaN, "s", 
            "The viscoelastic time step in seconds for the second velocity step in seconds"
            ),
        timestep = ParameterFloat(
            10000.0, "s", 
            "Main model time step in seconds used for forecasting viscoelastic stress, solving heat "
            *"conduction equation, advecting markers and advancing model time. The main model time step "
            *"may be reduced after solving the Stokes-continuity equations to ensure that markers do not "
            *"advect beyond a user defined limit. The model time step is reset to the viscoelastic time "
            *"step at the start of each time loop iteration. The heat equation loop may use time steps "
            *"that are smaller than the main model time step if temperature changes are too large"
            ),
        iupdate_timestep = ParameterInt(
            0, "None", 
            "0 = off; 1 = update timestep using maximum displacement"),

        # Timestep - Multiple increase parameters
        iuse_multiple_timestep_increase = ParameterInt(
            0, "None", 
            "Flag to use multiple viscoelastic time step increases: 0 = off; 1 = on. When active the "
            *"viscoelastic time step will be increased multiple times during the simulation at the time "
            *"steps specified by `ntime_increase_1`, `ntime_increase_2`, etc"),
        ntime_increase_1 = ParameterInt(
            100, "None", "First time step number at which to increase the timestep"),
        ntime_increase_2 = ParameterInt(
            100, "None", "Second time step number at which to increase the timestep"),
        ntime_increase_3 = ParameterInt(
            100, "None", "Third time step number at which to increase the timestep"),
        ntime_increase_4 = ParameterInt(
            100, "None", "Fourth time step number at which to increase the timestep"),
        time_increase_factor = ParameterFloat(
            2.0, "None", "Multiplication factor used to increase the "
            *"viscoelastic time step"
            ),
        cell_displ_factor = ParameterFloat(
            2.0, "None", 
            "Factor used to adjust the maximum allowed cell displacement fraction during multiple "
            *"viscoelastic time step increases"
            ),

        # Timestep - Thermal loop parameters
        timestep_heat = ParameterFloat(
            0.0, "s", "Thermal time step in seconds used in adaptive heat solver loop"),
        timestep_sum = ParameterFloat(
            0.0, "s", "Total thermal loop time in seconds"),

        # Timestep - Output steps parameters
        timestep_out = ParameterFloat(
            1.0e3, "s", "Time step for output"),
        nskip = ParameterInt(
            0, "None", 
            "Number of time steps to skip before output is generated. This is calculated as "
            *"`floor(Int, timestep_out/timestep_viscoelastic)` and is only used if `iuse_fixed_output_counter` is 1."
            ),
        icount_output = ParameterInt(
            0, "None", 
            "Output loop counter used to track time steps between outputs if `iuse_fixed_output_counter` is 1, "
            *"An output event is triggered when this counter is equal to nskip, which involves setting the output counter to zero. "
            *"This counter is incremented by one at the end of each time step."
            ),
        noutput = ParameterInt(
            0, "None", "Number of outputs generated during a model run"),
        time_of_next_output_myr = ParameterFloat(
            0.0, "Myr", "Time of next output in Myr"),
        # Timestep - Main time loop parameters
        iuse_fixed_output_counter = ParameterInt(
            1, "None",
            "Flag to use fixed output counter: 0 = off; 1 = on. When active, the model will output "
            *"at a fixed number of time steps based on an initial estimate at the start of the simulation. "
            *"If set to 0, then the model will output once total model time (timesum) reaches the next output time based "
            *"on the output time step. "
            *"Note: this flag is automatically forced to 0 during time-loop initialization when "
            *"`iuse_local_adaptive_time_stepping = 1`, since under local adaptive time stepping the "
            *"time step can vary significantly and the fixed-step counter (`nskip`) computed from the "
            *"initial time step is no longer a reliable output trigger."
            ),

        )
end
