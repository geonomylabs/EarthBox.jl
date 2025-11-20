"""
    ModelManager

Module responsible for managing the lifecycle and execution of the numerical model.

This module handles:
- Model initialization and setup
- Backup and restart functionality
- Model state management
- Time step progression
- Output management through XDMF files

The ModelManager acts as the central coordinator for the simulation, orchestrating
the interaction between various components like the grid, markers, solvers,
and output handlers.
"""
module ModelManager

include("backup/BackupManager.jl")
include("info/ModelInfo.jl")
include("xdmf/XdmfTimeStepsManager.jl")

import MPI
import EarthBox.PrintFuncs: @timeit_memit, print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.LoopInputStruct: LoopInput
import EarthBox: ConfigurationManager
import EarthBox.ConfigurationManager: Configuration
import EarthBox.EarthBoxDtypes: ParametersDictType
import EarthBox.ModelDataContainer.OutputStandard: OutputLists
import EarthBox: ModelParametersInfo
import EarthBox: StaggeredGrid
import EarthBox: MaterialGeometry
import EarthBox: Markers
import EarthBox: RockProperties
import EarthBox: Rheology
import EarthBox: Interpolation
import EarthBox: BoundaryConditions
import EarthBox: TimeStep
import EarthBox: StokesContinuitySolver
import EarthBox: HeatSolver
import EarthBox: GlobalPlasticityLoop
import EarthBox: Advection
import EarthBox: SurfaceProcesses
import EarthBox: MeltModel
import EarthBox: Serpentinization
import EarthBox: HydrothermalCirculation
import .XdmfTimeStepsManager: XdmfTimeSteps

"""
    ModelManagerState

State object for the model manager.

# Fields
- `restart_from_backup::Bool`
    - Controls whether to restart from a backup file called model_backup.jld 
       located in the model output directory.
- `paths::Dict{String, String}`
    - Dictionary of paths including model input file path (`model_input_file`), materials 
       input file path (`materials_input_file`), materials library file path 
       (`materials_library_file`), output directory path (`output_dir`), and storage 
       directory path (`storage_dir`).
- `run_paths::Dict{String, String}`
    - Dictionary of paths used to run the model including the model input file path 
       (`model_file`), materials input file path (`materials_file`), materials library 
       file path (`materials_library_file`), model output directory path (`output_dir`), 
       source directory (`src_dir`), and model output directory (`plots_dir`).
- `config::Configuration`
    - Configuration object for the model.
- `xdmf_time_steps::XdmfTimeSteps`
    - XDMF time steps object used to manage the export of XDMF files.
- `model::ModelData`
    - Model data object containing collections arrays and parameters.
- `output_lists::Union{OutputLists, Nothing}`
    - Output lists object used to manage the export of output files.
- `make_backup::Bool`
    - Controls whether to make a backup file called model_backup.jld 
       located in the model output directory after each output time step..
- `loop_time_start::Float64`
    - Start time of the loop.

# Constructor

    ModelManagerState(
        restart_from_backup::Bool,
        paths::Dict{String, String},
        run_paths::Dict{String, String},
        initialization_params::Union{ParametersDictType, Nothing},
        use_mumps::Bool,
        nprocs::Int,
        use_internal_mumps::Bool,
        mpi_comm::Union{MPI.Comm, Nothing},
        mpi_initialized::Bool,
        mpi_rank::Int,
        pass_large_arrays_via_mpi::Bool
    )

Main state constructor for the model manager.

# Inputs
- `restart_from_backup::Bool`
    - Controls whether to restart from a backup file called model_backup.jld 
       located in the model output directory.
- `paths::Dict{String, String}`
    - Dictionary of paths.
- `run_paths::Dict{String, String}`
    - Dictionary of paths used to run the model.
- `initialization_params::Union{ParametersDictType, Nothing}`
    - Dictionary of initialization parameters provided as inputs.
- `use_mumps::Bool`
    - Controls whether to use the MUMPS solver.
- `nprocs::Int`
    - Number of MPI processes to use for the MUMPS solver.
- `use_internal_mumps::Bool`
    - Controls whether to use the internal MUMPS solver..
- `mpi_comm::Union{MPI.Comm, Nothing}`
    - MPI communicator.
- `mpi_initialized::Bool`
    - Indicates whether MPI is initialized.
- `mpi_rank::Int`
    - The MPI rank of the current process.
- `pass_large_arrays_via_mpi::Bool`
    - Controls whether to pass large arrays via MPI between a parent and spawned 
       child process. If false arrays are passed via file IO. Passing large arrays 
       via MPI between a parent and spawned child process does not work on all systems 
       and is therefore tuned off by default.
"""
mutable struct ModelManagerState
    restart_from_backup::Bool
    paths::Dict{String, String}
    run_paths::Dict{String, String}
    config::Configuration
    xdmf_time_steps::XdmfTimeSteps
    model::ModelData
    output_lists::Union{OutputLists, Nothing}
    make_backup::Bool
    loop_time_start::Float64
end

function ModelManagerState(
    restart_from_backup::Bool,
    paths::Dict{String, String},
    run_paths::Dict{String, String},
    initialization_params::Union{ParametersDictType, Nothing},
    use_mumps::Bool,
    nprocs::Int,
    use_internal_mumps::Bool,
    mpi_comm::Union{MPI.Comm, Nothing},
    mpi_initialized::Bool,
    mpi_rank::Int,
    pass_large_arrays_via_mpi::Bool
)
    model = @timeit_memit "Finished initializing model" begin 
        ModelData(paths["model_input_file"], initialization_params)
    end
    
    return ModelManagerState(
        restart_from_backup,
        paths,
        run_paths,
        Configuration(
            run_paths, use_mumps, nprocs, use_internal_mumps, mpi_comm, 
            mpi_initialized, mpi_rank, pass_large_arrays_via_mpi
            ),
        XdmfTimeSteps(paths["output_dir"]),
        model,
        nothing,
        false,
        0.0
    )
end

function initialize_model!(manager::ModelManagerState)
    StaggeredGrid.initialize!(manager.model)
    initialize_geometry(manager.model)
    initialize_and_reset_boundary_conditions!(manager.model)
    Markers.initialize!(manager.model, manager.paths)
    initialize_rock_properties(manager.model)
    StokesContinuitySolver.initialize!(manager.model)
    Rheology.initialize!(manager.model)
    GlobalPlasticityLoop.initialize!(manager.model)
    HeatSolver.initialize!(manager.model)
    Advection.initialize!(manager.model)
    Interpolation.initialize!(manager.model)
    initialize_melt_model(manager.model)
    initialize_surface_processes(manager.model)
    TimeStep.MultipleTimeStepIncrease.initialize!(manager.model)
    TimeStep.SingleTimeStepIncrease.initialize!(manager.model)
    check_model_initialization(manager.model)
end

function initialize_geometry(model::ModelData)
    # EarthLayering is initialized to calculate derived parameters for
    # layer thicknesses
    MaterialGeometry.EarthLayering.initialize!(model)
    # LithoStrongZones is initialized to calculate derived parameters for 
    # inner boundaries of lithosphere strong zones based on T-type refinement
    MaterialGeometry.LithoStrongZones.initialize!(model)
end

function initialize_and_reset_boundary_conditions!(model::ModelData)
    BoundaryConditions.initialize!(model)
    BoundaryConditions.Temperature.initialize!(model)
    # TransientBottomTemperature is initialized to calculate derived parameters for
    # transient bottom temperature boundary condition based on provided inputs
    BoundaryConditions.TransientBottomTemperature.initialize!(model)
    BoundaryConditions.Velocity.initialize!(model)
    BoundaryConditions.Pressure.initialize!(model)
    BoundaryConditions.VelocityStep.initialize!(model)
    BoundaryConditions.VelocityStop.initialize!(model)
    BoundaryConditions.VelocityFromStrainRate.initialize!(model)
    # The reset bc function is not called when initializing the model via the API
    # which leads to vyu being inconsistent in parameters.yml files
    # This doesn't affect the model run since the boundary conditions are recalculated
    # at the start of each time step.
    BoundaryConditions.reset_bc!(model)
end

function initialize_rock_properties(model::ModelData)
    RockProperties.initialize!(model)
    HydrothermalCirculation.initialize!(model)
    Serpentinization.initialize!(model)
end
 
function initialize_melt_model(model::ModelData)
    MeltModel.initialize!(model)
    MeltModel.MeltDamage.initialize!(model)
    MeltModel.Extraction.initialize!(model)
    MeltModel.Extrusion.initialize!(model)
end

function initialize_surface_processes(model::ModelData)
    SurfaceProcesses.Topography.initialize!(model)
    SurfaceProcesses.Sealevel.initialize!(model)
    SurfaceProcesses.SedimentTransport.initialize!(model)
    SurfaceProcesses.SaltDeposition.initialize!(model)
end
 
function check_model_initialization(model::ModelData)
    SurfaceProcesses.Topography.check_topography_model(model)
end
 
function initialize_model_run(manager::ModelManagerState)
    if manager.restart_from_backup
        restart(manager)
    end
    set_output_lists(manager)
end
 
function restart(manager::ModelManagerState)
    BackupManager.LoadBackup.load_backup_jld2(
        manager.model.obj_dict, manager.paths["output_dir"])
end
 
function set_output_lists(manager::ModelManagerState)
    manager.output_lists = OutputLists(
        manager.model.markers,
        manager.model.heat_equation,
        manager.model.stokes_continuity
    )
end
 
function print_model_options(manager::ModelManagerState)
    model = manager.model
    StaggeredGrid.print_option(model)
    BoundaryConditions.print_option(model)
    Markers.MarkerMaterials.print_option(model)
    Markers.MarkerTemperature.print_option(model)
    RockProperties.print_option(model)
    StokesContinuitySolver.VelocityType.print_option(model)
    Rheology.PlasticFailure.MarkerPlasticity.print_option(model)
    GlobalPlasticityLoop.print_option(model)
    Advection.print_option(model)
end
 
function print_model_info(manager::ModelManagerState)
    print_number_of_markers(manager)
    print_number_of_stokes_unknowns(manager)
    print_number_of_heat_unknowns(manager)
    print_number_of_time_steps(manager)
    print_material_id_dicts(manager)
    # ModelParametersInfo.print_parameters(manager.model)
    print_model_options(manager)
    Rheology.FlowViscosity.print_model_flow_type_information(manager.model)
end
 
function print_number_of_markers(manager::ModelManagerState)
    marknum = manager.model.markers.parameters.distribution.marknum.value
    print_info("Number of markers: $(marknum)")
end

function print_number_of_stokes_unknowns(manager::ModelManagerState)
    nstokes = manager.model.stokes_continuity.parameters.build.N.value
    print_info("Number of Stokes unknowns: $(nstokes)")
end

function print_number_of_heat_unknowns(manager::ModelManagerState)
    nheat = manager.model.heat_equation.parameters.build.N.value
    print_info("Number of heat unknowns: $(nheat)")
end

function print_number_of_time_steps(manager::ModelManagerState)
    ntimestep_max = manager.model.timestep.parameters.main_time_loop.ntimestep_max.value
    print_info("Maximum number of time steps: $(ntimestep_max)")
end

function print_material_id_dicts(manager::ModelManagerState)
    mats = manager.model.materials.dicts
    print_info("", level=1)
    print_info("Material IDs for Material Domains")
    print_info("------------------------------------", level=1)
    for (key, value) in mats.matid_domains
        print_info("$(key) $(value)", level=2)
    end
    print_info("", level=1)
    print_info("Material IDs for Material Types", level=1)
    print_info("-------------------------------", level=1)
    for (key, value) in mats.matid_types
        print_info("$(key) $(value)", level=2)
    end
end

function export_array_and_parameter_info(manager::ModelManagerState)
    ModelInfo.export_collection_information(
        manager.model, manager.run_paths["output_dir"], "arrays")
    ModelInfo.export_collection_information(
        manager.model, manager.run_paths["output_dir"], "parameters")
end

function manage_model_output(manager::ModelManagerState)
    icount_output = manager.model.timestep.parameters.output_steps.icount_output.value
    make_files = manager.config.output.general["make_files"]
    if icount_output == 0 && make_files
        if manager.make_backup
            BackupManager.MakeBackup.make_backup_jld2(
                manager.model.obj_dict, manager.paths["output_dir"])
        end
        XdmfTimeStepsManager.export_xdmf(
            manager.xdmf_time_steps, manager.model, 
            manager.output_lists, manager.config.output
        )
    end
end

function get_loop_input_struct(
    model_manager::ModelManagerState
)::LoopInput
    bools = BoundaryConditions.get_boolean_options(model_manager.model)
    no_yielding_in_mobile_wall = get(bools, :no_yielding_in_mobile_wall, false)
    no_yielding_in_plate_extension = get(bools, :no_yielding_in_plate_extension, false)
    return LoopInput(
        model_manager.run_paths,
        model_manager.config.solver,
        no_yielding_in_mobile_wall,
        no_yielding_in_plate_extension
    )
end

function print_output_config(manager::ModelManagerState)
    ConfigurationManager.OutputConfig.print_output_config(manager.config.output)
end

end # module