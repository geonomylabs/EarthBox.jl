""" Run model

    run_time_steps(earthbox::EarthBoxState; kwargs...)

Initialize time step parameters and run the model for a given number of time steps or model duration.

# Arguments
- `earthbox::EarthBoxState`: EarthBoxState object

# Keyword arguments
- `make_backup::Bool`
    - Whether to make a backup of the model state at each output time step.
- `$(PDATA.model_duration_myr.name)::Float64`
    - $(PDATA.model_duration_myr.description)
- `$(PDATA.ntimestep_max.name)::Int`
    - $(PDATA.ntimestep_max.description)
- `$(PDATA.timestep_viscoelastic.name)::Float64`
    - $(PDATA.timestep_viscoelastic.description)
- `$(PDATA.timestep_out.name)::Float64`
    - $(PDATA.timestep_out.description)
- `$(PDATA.displ_limit.name)::Float64`
    - $(PDATA.displ_limit.description)
- `$(PDATA.iuse_boundary_displacement.name)::Int`
    - $(PDATA.iuse_boundary_displacement.description)
- `$(PDATA.strain_limit.name)::Float64`
    - $(PDATA.strain_limit.description)
- `$(PDATA.iuse_extensional_strain.name)::Int`
    - $(PDATA.iuse_extensional_strain.description)
- `$(PDATA.iupdate_timestep.name)::Int`
    - $(PDATA.iupdate_timestep.description)
- `$(PDATA.number_of_transport_timesteps_per_model_timestep.name)::Int`
    - $(PDATA.number_of_transport_timesteps_per_model_timestep.description)
- `$(PDATA.iuse_fixed_output_counter.name)::Int`
    - $(PDATA.iuse_fixed_output_counter.description)
"""
function run_time_steps(
    earthbox::EarthBoxState; 
    kwargs...
)::Nothing
    mpi_rank = earthbox.model_manager.config.solver.mpi_rank
    if mpi_rank == 0
        print_info("Number of Julia threads: $(Threads.nthreads())")
        print_info("Running time steps on mpi rank $mpi_rank...")
        make_backup = get(kwargs, :make_backup, false)
        earthbox.model_manager.make_backup = make_backup
        
        ModelManager.print_output_config(earthbox.model_manager)
        ModelManager.initialize_model_run(earthbox.model_manager)
        if earthbox.restart_from_backup
            # Do not use all inputs for time loop initialization if restarting
            # from backup since this will override timestep parameters
            # loaded from backup file. Only initialize number of time steps
            # and make backup since the user may want to change these
            # parameters when restarting from backup.
            TimeLoop.initialize!(
                earthbox.model_manager.model,
                ntimestep_max=get(kwargs, :ntimestep_max, nothing)
                )
        else
            TimeLoop.initialize!(
                earthbox.model_manager.model,
                model_duration_myr=get(kwargs, :model_duration_myr, nothing),
                ntimestep_max=get(kwargs, :ntimestep_max, nothing),
                timestep_viscoelastic=get(kwargs, :timestep_viscoelastic, nothing),
                timestep_out=get(kwargs, :timestep_out, nothing),
                iupdate_timestep=get(kwargs, :iupdate_timestep, nothing),
                iuse_boundary_displacement=get(kwargs, :iuse_boundary_displacement, nothing),
                iuse_extensional_strain=get(kwargs, :iuse_extensional_strain, nothing),
                displ_limit=get(kwargs, :displ_limit, nothing),
                strain_limit=get(kwargs, :strain_limit, nothing),
                number_of_transport_timesteps_per_model_timestep=get(
                    kwargs, :number_of_transport_timesteps_per_model_timestep, 5),
                iuse_fixed_output_counter=get(kwargs, :iuse_fixed_output_counter, nothing)
            )
        end
        ModelManager.print_model_info(earthbox.model_manager)
        ModelManager.export_array_and_parameter_info(earthbox.model_manager)
        paths_keys = earthbox.paths.key_names

        file_mover = FileMover.FileMoverState(
            output_directory=earthbox.paths.paths[paths_keys.output_dir],
            storage_directory=earthbox.storage_dir
        )
        try
            TimeLoop.run_loop!(earthbox.model_manager)
        catch e
            println("!!! ERROR !!! Type: $(typeof(e)) with message: $(e)")
            showerror(stdout, e, catch_backtrace())
            RunMumpsSolverLoop.shutdown_persistent_solver(
                earthbox.model_manager.config.solver)
            if MPIManager.is_mpi_initialized()
                MPIManager.finalize_mpi!()
            end
        finally
            FileMover.file_mover_cleanup(file_mover)
            if MPIManager.is_mpi_initialized()
                RunMumpsSolverLoop.shutdown_persistent_solver(
                    earthbox.model_manager.config.solver)
                MPIManager.finalize_mpi!()
            end
        end
    end

    return nothing
end

""" 
    run_model(; kwargs...)

Run EarthBox model using user specified paths.

# Keyword arguments
- `make_backup::Bool`
    - Whether to make a backup of the model state at each output time step.
- `restart_from_backup::Bool`
    - Whether to restart from a backup file.
- `use_mumps::Bool`
    - Whether to use MUMPS for solving the linear system of equations.
- `use_mumps_internal::Bool`
    - Whether to use the internal MUMPS solver.
- `nprocs::Int`
    - Number of MPI processes to use.
- `marker_output_from_user::Union{Dict{String, Bool}, Nothing}`
    - Dictionary of marker output keys and values from user input. See list of 
       marker output keys below.
- `model_input_file::String`
    - Path to the model input file.
- `materials_input_file::String`
    - Path to the materials input file.
- `materials_library_file::String`
    - Path to the materials library file.
- `model_output_path::Union{String, Nothing}`
    - Path to the model output directory. If not provided, the model output 
       directory is set to "output" and is created in the current working 
       directory.

# List of marker output keys
$(join(["- $(getfield(MarkerOutputKeys(), f))" for f in fieldnames(typeof(MarkerOutputKeys()))], "\n"))

# Returns
- `nothing`

"""
function run_model(;
    make_backup::Bool         = false,
    restart_from_backup::Bool = false,
    use_mumps::Bool           = false,
    use_mumps_internal        = true,
    nprocs::Int                = 1,
    marker_output_from_user::Union{Dict{String, Bool}, Nothing} = nothing,
    model_input_file::String,
    materials_input_file::String,
    materials_library_file::String,
    model_output_path::Union{String, Nothing} = nothing
)::Nothing
    if model_output_path === nothing
        model_output_path = "output"
    end
    println("model_input_file: $(model_input_file)")
    println("materials_input_file: $(materials_input_file)")
    println("materials_library_file: $(materials_library_file)")
    println("model_output_path: $(model_output_path)")
    check_cl_input_files(
        [model_input_file, materials_input_file, materials_library_file])
    eb = EarthBoxState(
        restart_from_backup = restart_from_backup,
        use_mumps           = use_mumps,
        nprocs              = nprocs,
        use_internal_mumps  = use_mumps_internal,
        paths               = Dict(
            "model_input_file"       => model_input_file,
            "materials_input_file"   => materials_input_file,
            "materials_library_file" => materials_library_file,
            "output_dir"             => model_output_path
        )
    )
    if marker_output_from_user !== nothing
        set_marker_output_keys!(eb, marker_output_from_user)
    end
    PRINT_SETTINGS.print_performance = false
    eb.model_manager.config.solver.verbose_output = 0
    ModelManager.initialize_model!(eb.model_manager)
    run_time_steps(eb, make_backup=make_backup)
    return nothing
end

"""
    set_marker_output_keys!(
        eb::EarthBoxState, 
        marker_output_from_user::Dict{String, Bool}
    )::Nothing

Set marker output keys from user input.

# Arguments
- `eb::EarthBoxState`
    - EarthBoxState object
- `marker_output_from_user::Dict{String, Bool}`
    - Dictionary of marker output keys and values from user input

# Returns
- `nothing`

"""
function set_marker_output_keys!(
    eb::EarthBoxState,
    marker_output_from_user::Dict{String, Bool}
)::Nothing
    marker_output = eb.model_manager.config.output.marker_output
    for (key, value) in marker_output_from_user
        if haskey(marker_output, key)
            marker_output[key] = value
            print_info("Marker output key $(key) set to $(value)")
        end
    end

    return nothing
end

function check_cl_input_files(
    input_files::Vector{String}
)::Nothing
    for input_file in input_files
        if !isfile(input_file)
            throw(ArgumentError("Input file not found: $(input_file)"))
        end
        if !endswith(input_file, ".yml")
            throw(ArgumentError("Input file must have .yml extension: $(input_file)"))
        end
    end
    return nothing
end

function greet()
    println("Hello, EarthBox!")
end
