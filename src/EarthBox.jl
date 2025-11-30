"""
    EarthBox

Main module for the EarthBox project.

"""
module EarthBox

# Load startup configuration
include("__init__.jl")

import MPI
import Base.catch_backtrace
import .ModelManager
import .StaggeredGrid
import .TtypeCalculator
import .TimeLoop
import .FileMover
import .MPIManager
import .ParallelSolver: RunMumpsSolverLoop
import .UnitConversion: UnitConversionData
import .ConfigurationManager.OutputConfig: MarkerOutputKeys
# Data types
import .EarthBoxDtypes: ParametersDictType, MaterialsDictType, MaterialDictType
import .ParameterRegistry: get_eb_parameters
# Utilities
import .PrintFuncs: print_info
import .PrintFuncs: print_warning, print_info, @timeit_memit, PRINT_SETTINGS
import .SysTools: get_username
import .GetModels: get_models
# Testing and benchmark management
import .TestManager
import .BenchmarksManager
# Material library management
import .MaterialLibraryCollection: MaterialLibrary
import .MaterialLibraryCollection: make_material_collection_table
import .MaterialLibraryCollection: get_material_names_list_string, make_material_table
import .MaterialLibraryCollection: get_material_collection
import .Markers.MarkerMaterials.Registry: MaterialTypesRegistry
import .Markers.MarkerMaterials.Registry: MaterialDomainsRegistry
# Case input management
import .CaseInputTools: get_case_inputs_and_convert_to_standard_units, initialize_cases,
    convert_case_parameters_to_standard_units!
import .CaseInputTools: CaseBuilder, CasePrinters, GetInputs
import .CaseInputTools.CaseTypes: CaseType, CaseCollectionType, CaseParameter
import .CaseInputTools.CaseBuilder: define_case_group!
import .CaseInputTools.CasePrinters: print_case_info
# Command-line argument management
import .GetArgs
import .GetArgs: get_case_name_from_cl_args, get_istart, get_iend, get_plot_option_name
import .GetArgs: get_model_output_path_from_args, get_earthbox_project_path_from_args
import .GetArgs: get_model_output_path_from_cl_args, get_runit_actions_from_cl_args
# Path management
import .EarthBoxPaths
import .GetPaths
import .GetPaths: get_model_output_path, get_storage_path
# Plotting management
import .PlotToolsManager: ModelPlots2DManager
import .PlotToolsManager.ModelPlots2DManager: ModelPlots2D, GridPlotsManager, 
    MarkerPlotsManager, PlasticityPlotManager, RheologyPlotsManager,
    plot_scalars, plot_markers, plot_velocity, plot_yield_strength, 
    plot_stokes_convergence, run_cl_plotter
import .PlotToolsManager.RunParallelPlot: run_parallel_marker_plotter
# RunTools management
import .RunTools: run_earthbox, remote_model_loop, local_model_loop
import .RunTools: execute_remote_script_in_background, execute_earthbox_script
# API Option name import
import .BenchmarksManager.option_names as benchmark_names
import .StaggeredGrid.option_names as grid_type_names
import .BoundaryConditions.option_names as bc_model_type_names
import .Markers.MarkerMaterials.option_names as material_model_names
import .Markers.MarkerTemperature.option_names as initial_temperature_names
import .Markers.MarkerCoordinates.option_names as marker_distribution_names
import .Markers.MarkerFriction.option_names as friction_init_names
import .RockProperties.RhoCpModel: option_names as rhocp_model_names
import .RockProperties.DensityModel: option_names as density_model_names
import .RockProperties.ThermalConductivityModel: option_names as thermal_conductivity_model_names
import .StokesContinuitySolver: velocity_type_names
import .Advection.option_names as advection_scheme_names
import .SurfaceProcesses.Topography.option_names as topo_node_advection_names
import .SurfaceProcesses.Sealevel.option_names as sealevel_option_names
import .GlobalPlasticityLoop.option_names as global_plasticity_names
import .SolidusLiquidus.SolidusModels: option_names as solidus_names
import .SolidusLiquidus.LiquidusModels: option_names as liquidus_names

# Test Export
export TestManager

# Option name export
export bc_model_type_names, velocity_type_names, advection_scheme_names, 
    topo_node_advection_names, sealevel_option_names, global_plasticity_names, 
    solidus_names, liquidus_names, grid_type_names, benchmark_names
## Rock property option names
export rhocp_model_names, density_model_names, thermal_conductivity_model_names
## Marker initialization option names
export material_model_names, marker_distribution_names, initial_temperature_names,
    friction_init_names

# Export API functions
export run_time_steps, print_info, run_model, get_models
# Export Data Structures and Type manager
export EarthBoxState, ModelDataContainer, get_eb_parameters, EarthBoxDtypes
# Export model management module
export ModelManager
# Export EarthBox API Modules
export StaggeredGrid, MaterialGeometry, BoundaryConditions, Markers, RockProperties
export HydrothermalCirculation, Serpentinization
export StokesContinuitySolver, GlobalPlasticityLoop, HeatSolver, Advection
export Interpolation, TimeStep, TimeLoop
export MeltModel, SurfaceProcesses, Gravity, BenchmarksManager
export Kinematics, Rheology, Regrid
# Export material library, types and domains modules 
export MaterialLibraryCollection, MaterialLibrary, MaterialTypesRegistry,
    MaterialDomainsRegistry, make_material_collection_table, 
    get_material_names_list_string, make_material_table, get_material_collection
export MaterialsDictType, MaterialDictType
# Export utility modules and functions
export ConversionFuncs, GetPaths, get_username
# Export PrintFuncs tools
export print_warning, print_info, @timeit_memit, PRINT_SETTINGS
# Export case management tools
export GetArgs
export get_case_inputs_and_convert_to_standard_units, CaseType 
export CaseCollectionType, CaseParameter
export convert_case_parameters_to_standard_units!
export GetInputs, print_case_info, initialize_cases, define_case_group!
# Export GetArgs tools
export get_case_name_from_cl_args, get_istart, get_iend, get_plot_option_name, 
    get_model_output_path_from_args, get_earthbox_project_path_from_args, 
    get_model_output_path_from_cl_args, get_runit_actions_from_cl_args
# Path tools
export get_model_output_path, get_storage_path
# Export plotting managers and functions
export ModelPlots2DManager, ModelPlots2D, GridPlotsManager, MarkerPlotsManager 
export PlasticityPlotManager, RheologyPlotsManager, run_parallel_marker_plotter 
export plot_scalars, plot_markers, plot_velocity, plot_yield_strength
export plot_stokes_convergence, run_cl_plotter
# Export RunTools tools
export run_earthbox, remote_model_loop, local_model_loop
export execute_remote_script_in_background, execute_earthbox_script
# Export marker output configuration keys
export MarkerOutputKeys
# Export GetModels tools
export get_models

const PDATA = get_eb_parameters()

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
       and "ParMETIS". The "ParMETIS" option. The "PTSCOTCH" is sometimes not
       available with Julia MUMPS binaries.
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
            analysis_method, parallel_ordering_method, memory_relax_perc, verbose_output
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

""" Make initialization parameter dictionary

The initialization_params dictionary provides the user a way to define key
parameters for initializing the model via the API including model dimensions,
nodal resolution and grid refinement parameters. If this information is not 
provided (i.e. initialization_params = nothing), the model will be initialized 
with parameters read from the model input file.

"""
function make_initialization_parameter_dict(
    xnum::Union{Int, Nothing},
    ynum::Union{Int, Nothing},
    znum::Union{Int, Nothing},
    dx_marker::Union{Float64, Nothing},
    dy_marker::Union{Float64, Nothing},
    dz_marker::Union{Float64, Nothing},
    nmarkers_cell_x::Union{Float64, Nothing},
    nmarkers_cell_y::Union{Float64, Nothing},
    nmarkers_cell_z::Union{Float64, Nothing},
    xsize::Union{Float64, Nothing},
    ysize::Union{Float64, Nothing},
    zsize::Union{Float64, Nothing},
    type_refinement_parameters::Union{Dict{String, Float64}, Nothing}
)::Union{ParametersDictType, Nothing}
    if type_refinement_parameters !== nothing
        initialization_params = define_init_parameters_for_ttype_grid(
            xsize, ysize, type_refinement_parameters
        )
    else
        initialization_params = define_init_parameters_for_general_basic_grid(
            xsize, ysize, zsize, xnum, ynum, znum
        )
    end
    add_marker_resolution_parameters!(
        initialization_params,
        dx_marker, dy_marker, dz_marker,
        nmarkers_cell_x, nmarkers_cell_y, nmarkers_cell_z
    )
    return initialization_params
end

function define_init_parameters_for_general_basic_grid(
    xsize::Union{Float64, Nothing},
    ysize::Union{Float64, Nothing},
    zsize::Union{Float64, Nothing},
    xnum::Union{Int, Nothing},
    ynum::Union{Int, Nothing},
    znum::Union{Int, Nothing}
)::Union{ParametersDictType, Nothing}
    # These minimum attributes are used to determine if the user has provided
    # the minimum number of input parameters to initialize the model. If the
    # has not provided the minimum number of input parameters nothing is returned.
    # This instructs the code to use parameters from a model input file as
    # opposed to using the user provided parameters via the API.
    minimum_attr_list = [xsize, ysize, xnum, ynum]
    none_count = count(==(nothing), minimum_attr_list)
    # If some but not all key initialization parameters are provided, throw an error as the user
    # has not provided sufficient inputs for initializing the model and ambiguity is introduced.
    if 0 < none_count < length(minimum_attr_list)
        throw(ArgumentError(
            "Sufficient API inputs were not provided for initializing the general basic grid:  "
            *"xsize = $xsize, ysize = $ysize, xnum = $xnum, ynum = $ynum. Either provide all "
            *"inputs via the API or define entries for these parameters in the model input file."
            )
            )
    end
    # If all key initialization parameters are provided (none_count = 0), return the 
    # initialization_params dictionary
    if none_count == 0
        print_info("Sufficient API inputs were provided for initializing the general basic grid.")
        initialization_params =  ParametersDictType(
            "xnum" => [Int(xnum), "None"],
            "ynum" => [Int(ynum), "None"],
            "xsize" => [Float64(xsize), "m"],
            "ysize" => [Float64(ysize), "m"]
        )
        # Add 3d parameters if they are provided by the user.
        if znum !== nothing
            initialization_params["znum"] = [Int(znum), "None"]
        end
        if zsize !== nothing
            initialization_params["zsize"] = [Float64(zsize), "m"]
        end
        return initialization_params
    else
        return nothing
    end
end

function define_init_parameters_for_ttype_grid(
    xsize::Union{Float64, Nothing},
    ysize::Union{Float64, Nothing},
    type_refinement_parameters::Dict{String, Float64}
)::Union{ParametersDictType, Nothing}
    (
        xo_highres, xf_highres, dx_highres, dx_lowres, 
        yf_highres, dy_highres, dy_lowres
    ) = TtypeCalculator.get_minimum_ttype_refinement_parameters(type_refinement_parameters)
    minimum_attr_list = [
        xsize, ysize, 
        xo_highres, xf_highres, dx_highres, dx_lowres, 
        yf_highres, dy_highres, dy_lowres
    ]
    # These minimum attributes are used to determine if the user has provided
    # the minimum number of input parameters to initialize the model. If the
    # has not provided the minimum number of input parameters nothing is returned.
    # This instructs the code to use parameters from a model input file as
    # opposed to using the user provided parameters via the API.
    none_count = count(==(nothing), minimum_attr_list)
    if 0 < none_count < length(minimum_attr_list)
        throw(ArgumentError(
            "Sufficient API inputs were not provided for initializing the T-type refined grid:  "
            *"xsize = $xsize, ysize = $ysize, xo_highres = $xo_highres, xf_highres = $xf_highres, "
            *"dx_highres = $dx_highres, dx_lowres = $dx_lowres, yf_highres = $yf_highres, "
            *"dy_highres = $dy_highres, dy_lowres = $dy_lowres. Either provide all "
            *"inputs via the API or define entries for these parameters in the model input file."
            )
            )
    end
    if none_count == 0
        print_info(
            "Parameters for initializing the T-type refined grid are provided as "
            * "input parameters from the API"
            )
        initialization_params =  ParametersDictType(
            "xsize" => [Float64(xsize), "m"],
            "ysize" => [Float64(ysize), "m"],
            "xo_highres" => [Float64(xo_highres), "m"],
            "xf_highres" => [Float64(xf_highres), "m"],
            "dx_highres" => [Float64(dx_highres), "m"],
            "dx_lowres" => [Float64(dx_lowres), "m"],
            "yf_highres" => [Float64(yf_highres), "m"],
            "dy_highres" => [Float64(dy_highres), "m"],
            "dy_lowres" => [Float64(dy_lowres), "m"]
        )
        return initialization_params
    else
        return nothing
    end
end

function add_marker_resolution_parameters!(
    initialization_params::Union{ParametersDictType, Nothing},
    dx_marker::Union{Float64, Nothing},
    dy_marker::Union{Float64, Nothing},
    dz_marker::Union{Float64, Nothing},
    nmarkers_cell_x::Union{Float64, Nothing},
    nmarkers_cell_y::Union{Float64, Nothing},
    nmarkers_cell_z::Union{Float64, Nothing}
)::Nothing

    if dx_marker !== nothing && initialization_params !== nothing
        print_info("Marker spacing in x direction is provided as input parameter from the API")
        initialization_params["dx_marker"] = [Float64(dx_marker), "m"]
    elseif nmarkers_cell_x !== nothing && initialization_params !== nothing
        print_info("Number of markers per cell in x direction is provided as input parameter from the API")
        initialization_params["nmarkers_cell_x"] = [Float64(nmarkers_cell_x), "None"]
    elseif nmarkers_cell_x == nothing && dx_marker == nothing && initialization_params !== nothing
        error(
            "You are trying to initialize the model through the API but have not provided either the "
            * "marker spacing in x direction (dx_marker) or the number of markers per cell in the "
            * "x-direction (nmarkers_cell_x). One of these parameters are required to initialize the "
            * "marker swarm."
            )
    end

    if dy_marker !== nothing && initialization_params !== nothing
        print_info("Marker spacing in y direction is provided as input parameter from the API")
        initialization_params["dy_marker"] = [Float64(dy_marker), "m"]
    elseif nmarkers_cell_y !== nothing && initialization_params !== nothing
        print_info("Number of markers per cell in y direction is provided as input parameter from the API")
        initialization_params["nmarkers_cell_y"] = [Float64(nmarkers_cell_y), "None"]
    elseif nmarkers_cell_y == nothing && dy_marker == nothing && initialization_params !== nothing
        error(
            "You are trying to initialize the model through the API but have not provided either the "
            * "marker spacing in y direction (dy_marker) or the number of markers per cell in the "
            * "y-direction (nmarkers_cell_y). One of these parameters are required to initialize the "
            * "marker swarm."
            )
    end
    # Add statement for 3D once implemented
    #if dz_marker !== nothing
    #    initialization_params["dz_marker"] = [Float64(dz_marker), "m"]
    #else
    #    initialization_params["nmarkers_cell_z"] = [Float64(nmarkers_cell_z), "None"]
    #end
    return nothing
end

function check_inputs(
    eb_paths::EarthBoxPaths.EarthBoxPathsState, 
    initialization_params::Union{ParametersDictType, Nothing}
)::Nothing
    if !haskey(eb_paths.paths, "model_input_file") && initialization_params === nothing
        throw(ArgumentError(
            "Both `model_input_file` and initialization parameters are not defined. EarthBox "
            *"cannot be initialized. One of these must be defined."
            ))
    end
    if haskey(eb_paths.paths, "model_input_file") && initialization_params !== nothing
        print_warning(
            "Model input file is provided and key initialization parameters are provided as input "
            *"parameters from the API. The model input file will be ignored."
            )
    end
    return nothing
end

function get_base_path()
    return @__DIR__
end

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
                    kwargs, :number_of_transport_timesteps_per_model_timestep, 5)
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

end # module