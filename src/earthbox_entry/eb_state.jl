"""
    EarthBoxState

Mutable struct used to store EarthBox model state.

# Fields
- `restart_from_backup::Bool`
    - Controls whether to restart model from a backup file called model_backup.jld 
       located in the output directory.
- `paths::Union{EarthBoxPaths.EarthBoxPathsState, Nothing}`
    - Struct containing path information including model input file, materials 
       input file, materials library file, output directory, and storage directory.
       See **Input Path Configuration** section below for more details.
- `xnum::Union{Int, Nothing}`
    - Number of basic nodes in the x direction
- `ynum::Union{Int, Nothing}`
    - Number of basic nodes in the y direction
- `xsize::Union{Float64, Nothing}`
    - Model width in meters
- `ysize::Union{Float64, Nothing}`
    - Model height in meters
- `nmarkers_cell_x::Union{Float64, Nothing}`
    - Number of markers per cell in the x direction
- `nmarkers_cell_y::Union{Float64, Nothing}`
    - Number of markers per cell in the y direction
- `initialization_params::Union{ParametersDictType, Nothing}`
    - Dictionary of initialization parameters provided as inputs
- `run_paths::Union{Dict{String, String}, Nothing}`
    - Dictionary containing key paths used for running EarthBox models including
       the model input file path (`model_file`), materials input file path 
       (`materials_file`), materials library file path (`materials_library_file`), 
       model output directory path (`output_dir`), source directory (`src_dir`), 
       and model output directory (`plots_dir`).
- `storage_dir::Union{String, Nothing}`
    - Path to storage directory where model backup files will be moved for 
       longer term storage
- `model_manager::Union{ModelManager.ModelManagerState, Nothing}`
    - Model manager state that stores key structs for managing the model run 
       including including the ModelData struct composed of collections of arrays 
       and parameters.
- `units::Union{UnitConversionData, Nothing}`
    - Unit conversion data

# Constructor

    EarthBoxState(;
        restart_from_backup::Bool = false, 
        paths::Union{Dict{String, String}, Nothing} = nothing, 
        xnum::Union{Int, Nothing} = nothing, 
        ynum::Union{Int, Nothing} = nothing, 
        znum::Union{Int, Nothing} = nothing, 
        xsize::Union{Float64, Nothing} = nothing, 
        ysize::Union{Float64, Nothing} = nothing, 
        zsize::Union{Float64, Nothing} = nothing, 
        nmarkers_cell_x::Union{Float64, Nothing} = nothing, 
        nmarkers_cell_y::Union{Float64, Nothing} = nothing, 
        nmarkers_cell_z::Union{Float64, Nothing} = nothing, 
        dx_marker::Union{Float64, Nothing} = nothing, 
        dy_marker::Union{Float64, Nothing} = nothing, 
        dz_marker::Union{Float64, Nothing} = nothing, 
        ttype_refinement_parameters::Union{Dict{String, Float64}, Nothing} = nothing, 
        use_mumps::Bool = false, 
        nprocs::Int = 1, 
        use_internal_mumps::Bool = true
        )

Main state constructor for EarthBox geodynamic simulations used to initialize 
the model and manage model run actions. The creation of this state object is 
accompanied by the initialization of model data structures including arrays and 
parameters. MPI is also initialized if the MUMPS solver is used.

# Inputs
- `restart_from_backup::Bool`
    - Controls whether to restart model from a backup file called model_backup.jld 
       located in the output directory defined in `paths`.
- `paths::Union{Dict{String, String}, Nothing}`
    - Dictionary of user provided input paths including model input file, materials 
       input file, materials library file, output directory, and storage directory. 
       See **Input Path Configuration** section below for more details.
- `xnum::Union{Int, Nothing}`
    - Number of basic nodes in the x direction
- `ynum::Union{Int, Nothing}`
    - Number of basic nodes in the y direction
- `znum::Union{Int, Nothing}`
    - Number of basic nodes in the z direction
- `xsize::Union{Float64, Nothing}`
    - Model width in meters
- `ysize::Union{Float64, Nothing}`
    - Model height in meters
- `zsize::Union{Float64, Nothing}`
    - Model depth in meters
- `nmarkers_cell_x::Union{Float64, Nothing}`
    - Number of markers per cell in the x direction
- `nmarkers_cell_y::Union{Float64, Nothing}`
    - Number of markers per cell in the y direction
- `nmarkers_cell_z::Union{Float64, Nothing}`
    - Number of markers per cell in the z direction
- `dx_marker::Union{Float64, Nothing}`
    - Marker spacing (meters) in the x direction used to calculate `nmarkers_cell_x`
- `dy_marker::Union{Float64, Nothing}`
    - Marker spacing (meters) in the y direction used to calculate `nmarkers_cell_y`
- `dz_marker::Union{Float64, Nothing}`
    - Marker spacing (meters) in the z direction used to calculate `nmarkers_cell_z`
- `ttype_refinement_parameters::Union{Dict{String, Float64}, Nothing}`
    - Dictionary of T-type refinement parameters. See the 
       **T-type Refinement Parameters** section below for more details.
- `use_mumps::Bool`
    - Controls whether to use the MUMPS solver.
- `nprocs::Int`
    - Number of MPI processes to use for the MUMPS solver.
- `analysis_method::Union{String, Symbol, Nothing}`
    - Analysis method to use for the MUMPS solver. Options are "PARALLEL" and "SERIAL".
       Default is "SERIAL". Using the "PARALLEL" option can significantly speed up the 
       MUMPS solver especially for large systems of equations. However, the "PARALLEL" 
       option may require MPI.jl to be configured using system libraries, MUMPS to be
       locally compiled with the appropriate system libraries and MUMPS.jl to be
       reconfigured. See installation instructions in the docs for more details.
       Set to "SERIAL" if you are having issues with the "PARALLEL" option.
- `parallel_ordering_method::Union{String, Symbol, Nothing}`
    - Parallel ordering method to use for the MUMPS solver. Options are "PTSCOTCH" 
       and "ParMETIS". The "ParMETIS" option is the default. The "PTSCOTCH" option 
       is sometimes not available with Julia MUMPS binaries.
- `memory_relax_perc::Int`
    - Memory relaxation percentage to use for the MUMPS solver. Default is 25.
- `verbose_output::Int`
    - Verbose output level to use for the MUMPS solver. Default is 0.
- `use_internal_mumps::Bool`
    - Controls whether to use the internal MUMPS solver. See the 
       **Internal vs External MUMPS Solver** section below for more details.
- `pass_large_arrays_via_mpi::Bool`
    - Controls whether to pass large arrays via MPI between a parent and spawned 
       child process. If false arrays are passed via file IO. Passing large arrays 
       via MPI between a parent and spawned child process does not work on all systems 
       and is therefore tuned off by default.

# Outputs
- `EarthBoxState` -- EarthBoxState object

# Input-dependent Initialization Scenarios

All model arrays are pre-allocated during initialization which requires `xnum`, 
`ynum`, `xsize`, and `ysize` to be defined for grid-related arrays and 
`nmarkers_cell_x` and `nmarkers_cell_y` to be defined for marker-related arrays. 
The procedure used to define these key initialization parameters is dependent
on the inputs provided by the user. For cases involving local grid refinement 
(e.g. T-type refinement), `xnum`, `ynum` are derived from the T-type refinement 
parameters. If marker spacing parameters are provided, `nmarkers_cell_x` and 
`nmarkers_cell_y` are derived from the marker spacing and `xsize` and `ysize` and 
corresponding average marker spacing parameters, `dx_marker` and `dy_marker`, 
respectively. Specific, input-dependent initialization scenarios are described below.

## Initialization Scenarios for Scalar Arrays

**Scenario 1**: Input parameters `xnum`, `ynum`, `xsize`, and `ysize` are defined 
as non-nothing API inputs:
- `xnum`, `ynum`, `xsize`, and `ysize` will be used to initialize the model 
    including scalar arrays

**Scenario 2**: `ttype_refinement_parameters`, `xsize`, and `ysize` are 
provided as non-nothing API inputs:
- `xnum` and `ynum` are automatically derived from the dictionary parameters and 
    used to initialize the model including scalar arrays

**Scenario 3**: A model input file containing entries for `xnum`, `ynum`, `xsize`, 
and `ysize` is defined in the `paths` dictionary and `ttype_refinement_parameters`, 
`xsize`, and `ysize` are equal to nothing:
- `xnum`, `ynum`, `xsize`, and `ysize` will be read from the model input file and used to initialize 
    the model including scalar arrays

Other scenarios (e.g. incomplete inputs) will throw an error. For example, if 
the user provides only `xnum`, and other parameters are equal to nothing, the 
code will throw an error as inputs are incomplete.

If a complete input scenario is provided and a path to a model input file is 
provided, the code will ignore the model input file defined in the `paths` 
dictionary and all inputs must be defined via the API with subsequent API 
initialization calls.

## Initialization Scenarios for Marker Arrays

**Scenario 1**: Input parameters `nmarkers_cell_x` and `nmarkers_cell_y` are all 
defined as inputs
- `nmarkers_cell_x` and `nmarkers_cell_y` will be used to initialize marker arrays
    
**Scenario 2**: Input parameters `dx_marker` and `dy_marker` are all defined as 
non-nothing inputs
- `nmarkers_cell_x` and `nmarkers_cell_y` will be derived from the marker spacing and `xsize` 
    and `ysize`
    
**Scenario 3**: A model input file is defined in the `paths` dictionary and 
`nmarkers_cell_x`, `nmarkers_cell_y` are defined in the model input file:
- `nmarkers_cell_x` and `nmarkers_cell_y` will be read from the model input file and used to 
    initialize marker arrays

**Scenario 4**: A model input file is defined in the `paths` dictionary and `dx_marker`, and 
`dy_marker` and `xsize` and `ysize` are defined in the model input file:
- `nmarkers_cell_x` and `nmarkers_cell_y` will be derived from the marker spacing and `xsize` 
    and `ysize`

If parameters `dx_marker` and `dy_marker` are encountered by EarthBox they will be 
used to calculate `nmarkers_cell_x` and `nmarkers_cell_y`. Marker resolution 
parameters in the x- and y-directions are calculated independently of each other.

# Input Path Configuration (paths dictionary)

The `paths` dictionary with type Dict{String, String} accepts the following optional key:value pairs:
- `model_input_file`
    - Path to yaml model input file containing input parameters for the model. 
       If omitted the user will have to define input parameters through various 
       initialization methods in the EarthBox API.

- `materials_input_file`
    - Path to yaml materials input file containing material information for each 
       material used in the model. If omitted the user will have to provide provide 
       this path or a dictionary of material inputs as inputs for marker material 
       initialization through the EarthBox API.

- `materials_library_file`
    - Path to materials library file containing material parameters associated with 
       the material names used in the materials input file. If omitted the user will 
       have to provide the path to the materials library file as input for marker 
        material initialization through the EarthBox API. A collection of material 
        libraries are provided in the EarthBox package [MaterialLibrary](@ref).

- `output_dir`
    - Path to output directory where model output files will be written. If omitted 
       the output directory will be set to the active model directory with the 
       default name 'output_dir'.

- `storage_dir`
    - Path to storage directory where output files will be moved for longer term 
       storage. If omitted output files will not be moved. omitted output files 
       will not be moved.

# T-type Refinement Parameters

The `ttype_refinement_parameters` dictionary with type Dict{String, Float64} 
accepts the following key:value pairs:
- `xo_highres`: Starting x-coordinate of high resolution region in meters
- `xf_highres`: Ending x-coordinate of high resolution region in meters
- `dx_highres`: Constant grid spacing in x-direction for high resolution region in meters
- `yf_highres`: Ending y-coordinate of high resolution region in meters
- `dy_highres`: Constant grid spacing in y-direction for high resolution region in meters
- `dx_lowres`: Average grid spacing in x-direction for low resolution regions in meters
- `dy_lowres`: Average grid spacing in y-direction for low resolution regions in meters

If all keys are defined in the `ttype_refinement_parameters` dictionary, the number of basic nodes in 
the x and y directions, `xnum` and `ynum`, will be calculated from the T-type refinement parameters.

Note that `dx_lowres` and `dy_lowres` are used to calculate the number of basic nodes in x and y directions
during model initialization and are not used in grid generation functions. Furthermore, parameters 
`dx_lowres` and `dy_lowres` are not required when initializing the staggered grid 
through [`StaggeredGrid.initialize!`](@ref ModelManager.StaggeredGrid.initialize!).

Parameters from the `ttype_refinement_parameters` dictionary are loaded into the model data structures
by this API call eliminating the need to provide them when initializing the staggered grid through
the EarthBox API [`StaggeredGrid.initialize!`](@ref ModelManager.StaggeredGrid.initialize!).

# Internal vs External MUMPS Solver

The bool `use_internal_mumps` controls whether the internal MUMPS solver is used where a persistent 
MPI process is spawned and information is passed between a parent and child Julia process via MPI. 
If `use_internal_mumps` is false, the MUMPS solver is ran as an external process and information 
is passed via file IO. The internal solver is preferred since it significantly minimizes Julia 
compilation and loading overhead and eliminates overhead associated with file IO.

"""
mutable struct EarthBoxState
    restart_from_backup::Bool
    paths::Union{EarthBoxPaths.EarthBoxPathsState, Nothing}
    xnum::Union{Int, Nothing}
    ynum::Union{Int, Nothing}
    xsize::Union{Float64, Nothing}
    ysize::Union{Float64, Nothing}
    nmarkers_cell_x::Union{Float64, Nothing}
    nmarkers_cell_y::Union{Float64, Nothing}
    initialization_params::Union{ParametersDictType, Nothing}
    run_paths::Union{Dict{String, String}, Nothing}
    storage_dir::Union{String, Nothing}
    model_manager::Union{ModelManager.ModelManagerState, Nothing}
    units::Union{UnitConversionData, Nothing}
end

function EarthBoxState(; 
    restart_from_backup::Bool=false, 
    paths::Union{Dict{String, String}, Nothing}=nothing,
    xnum::Union{Int, Nothing}=nothing,
    ynum::Union{Int, Nothing}=nothing,
    znum::Union{Int, Nothing}=nothing,
    xsize::Union{Float64, Nothing}=nothing,
    ysize::Union{Float64, Nothing}=nothing,
    zsize::Union{Float64, Nothing}=nothing,
    nmarkers_cell_x::Union{Float64, Nothing}=nothing,
    nmarkers_cell_y::Union{Float64, Nothing}=nothing,
    nmarkers_cell_z::Union{Float64, Nothing}=nothing,
    dx_marker::Union{Float64, Nothing}=nothing,
    dy_marker::Union{Float64, Nothing}=nothing,
    dz_marker::Union{Float64, Nothing}=nothing,
    ttype_refinement_parameters::Union{Dict{String, Float64}, Nothing}=nothing,
    use_mumps::Bool=false,
    nprocs::Int=1,
    analysis_method::Union{String, Symbol, Nothing}=nothing,
    parallel_ordering_method::Union{String, Symbol, Nothing}=nothing,
    memory_relax_perc::Union{Int, Nothing}=nothing,
    verbose_output::Union{Int, Nothing}=nothing,
    use_internal_mumps::Bool=true,
    pass_large_arrays_via_mpi::Bool=false
)
    mpi_comm = MPIManager.initialize_mpi!(use_mumps, use_internal_mumps)
    mpi_initialized = MPIManager.is_mpi_initialized()
    mpi_rank = MPIManager.get_mpi_rank(mpi_comm)
    mpi_implementation = MPIManager.get_mpi_implementation()
    print_info("MPI implementation is system?: $(mpi_implementation.is_system_mpi)", level=1)

    eb_paths = nothing
    initialization_params = nothing
    run_paths = nothing
    storage_dir = nothing
    model_manager = nothing
    units = nothing
    if mpi_rank == 0
        print_info("Initializing EarthBox on mpi rank $mpi_rank...", level=1)
        # Dictionary of paths attribute used to define input and output paths
        eb_paths = EarthBoxPaths.EarthBoxPathsState()
        EarthBoxPaths.set_paths!(eb_paths, paths)
        initialization_params = make_initialization_parameter_dict(
            xnum, ynum, znum, 
            dx_marker, dy_marker, dz_marker,
            nmarkers_cell_x, nmarkers_cell_y, nmarkers_cell_z,
            xsize, ysize, zsize,
            ttype_refinement_parameters
        )
        check_inputs(eb_paths, initialization_params)
        # The run_paths dictionary contains key paths used for running EarthBox models
        # and derived paths.
        run_paths = EarthBoxPaths.create_run_paths_dict(eb_paths, get_base_path())
        # Storage directory where output files will be moved in the background
        # Files will not be moved if this attribute is None.
        storage_dir = get(eb_paths.paths, "storage_dir", nothing)
        model_manager = ModelManager.ModelManagerState(
            restart_from_backup,
            Dict{String,String}(eb_paths.paths), 
            Dict{String,String}(run_paths),
            initialization_params,
            use_mumps,
            nprocs,
            use_internal_mumps,
            mpi_comm,
            mpi_initialized,
            mpi_rank,
            pass_large_arrays_via_mpi
            )
        set_mumps_solver_options(
            model_manager, 
            analysis_method, parallel_ordering_method, 
            memory_relax_perc, verbose_output
        )
        RunMumpsSolverLoop.initialize_persistent_solver(model_manager.config.solver, mpi_comm)
        EarthBoxPaths.check_output_dir(eb_paths)
        EarthBoxPaths.print_input_path_info(eb_paths)
        units = UnitConversionData()
    end
    
    return EarthBoxState(
        restart_from_backup,
        eb_paths,
        xnum,
        ynum,
        xsize,
        ysize,
        nmarkers_cell_x,
        nmarkers_cell_y,
        initialization_params,
        run_paths,
        storage_dir,
        model_manager,
        units
    )
end

function set_mumps_solver_options(
    model_manager::ModelManager.ModelManagerState, 
    analysis_method::Union{String, Symbol, Nothing}, 
    parallel_ordering_method::Union{String, Symbol, Nothing}, 
    memory_relax_perc::Union{Int, Nothing}, 
    verbose_output::Union{Int, Nothing}
)
    if analysis_method !== nothing && in(String(analysis_method), ["PARALLEL", "SERIAL"])
        model_manager.config.solver.analysis_method = String(analysis_method)
    end
    if parallel_ordering_method !== nothing && in(String(parallel_ordering_method), ["PTSCOTCH", "ParMETIS"])
        model_manager.config.solver.parallel_ordering_method = String(parallel_ordering_method)
    end
    if memory_relax_perc !== nothing && 0 <= memory_relax_perc <= 100
        model_manager.config.solver.memory_relax_perc = memory_relax_perc
    end
    if verbose_output !== nothing && 0 <= verbose_output <= 3
        model_manager.config.solver.verbose_output = verbose_output
    end
end
