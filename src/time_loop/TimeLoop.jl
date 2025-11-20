module TimeLoop

include("stopping/Stopping.jl")

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.PrintFuncs: @timeit_memit, print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ModelDataContainer: set_parameter!
import EarthBox.LoopInputStruct: LoopInput
import EarthBox: GridFuncs
import EarthBox: StaggeredGrid
import EarthBox: Regrid
import EarthBox: ModelManager
import EarthBox.ModelManager: ModelManagerState
import EarthBox: BoundaryConditions
import EarthBox.Markers: MarkerFriction
import EarthBox: RockProperties
import EarthBox: Interpolation
import EarthBox: Rheology
import EarthBox: GlobalPlasticityLoop
import EarthBox: StokesContinuitySolver
import EarthBox: HeatSolver
import EarthBox: Advection
import EarthBox: SurfaceProcesses
import EarthBox: TimeStep
import EarthBox: MeltModel
import EarthBox: HydrothermalCirculation
import EarthBox: Serpentinization
import EarthBox.Arrays: ArrayUtils
import EarthBox.TimeStep: TimeStepCalculator
import .Stopping

const PDATA = get_eb_parameters()

struct ValidInputNames
    model_duration_myr::Symbol
    ntimestep_max::Symbol
    timestep_viscoelastic::Symbol
    timestep_out::Symbol
    iupdate_timestep::Symbol
    iuse_boundary_displacement::Symbol
    iuse_extensional_strain::Symbol
    displ_limit::Symbol
    strain_limit::Symbol
    number_of_transport_timesteps_per_model_timestep::Symbol
end

""" Initialize time loop

# Arguments
- `model::ModelData`
    - Model data object containing the model parameters and arrays.

# Keyword Arguments
- `$(PDATA.model_duration_myr.name)::Union{Float64, Nothing}`
    - $(PDATA.model_duration_myr.description)
- `$(PDATA.ntimestep_max.name)::Union{Int, Nothing}`
    - $(PDATA.ntimestep_max.description)
- `$(PDATA.timestep_viscoelastic.name)::Union{Float64, Nothing}`
    - $(PDATA.timestep_viscoelastic.description)
- `$(PDATA.timestep_out.name)::Union{Float64, Nothing}`
    - $(PDATA.timestep_out.description)
- `$(PDATA.iupdate_timestep.name)::Union{Int, Nothing}`
    - $(PDATA.iupdate_timestep.description)
- `$(PDATA.iuse_boundary_displacement.name)::Union{Bool, Nothing}`
    - $(PDATA.iuse_boundary_displacement.description)
- `$(PDATA.displ_limit.name)::Union{Float64, Nothing}`
    - $(PDATA.displ_limit.description)
- `$(PDATA.iuse_extensional_strain.name)::Union{Bool, Nothing}`
    - $(PDATA.iuse_extensional_strain.description)
- `$(PDATA.strain_limit.name)::Union{Float64, Nothing}`
    - $(PDATA.strain_limit.description)
- `$(PDATA.number_of_transport_timesteps_per_model_timestep.name)::Union{Int, Nothing}`
    - $(PDATA.number_of_transport_timesteps_per_model_timestep.description)
"""
function initialize!(
    model::ModelData;
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    TimeStepCalculator.update_time_step_parameters!(model)
    TimeStepCalculator.update_transport_timestep_for_symmetric_extension!(
        model, get(kwargs, :number_of_transport_timesteps_per_model_timestep, 5))
    update_model_time_step!(model)
    update_nskip!(model)
    return nothing
end

function update_nskip!(model::ModelData)::Nothing
    timestep_params = model.timestep.parameters
    timestep_out = timestep_params.output_steps.timestep_out.value
    timestep_viscoelastic = timestep_params.main_time_loop.timestep_viscoelastic.value
    nskip = floor(Int, timestep_out/timestep_viscoelastic)
    set_parameter!(model, "nskip", nskip)
    return nothing
end

function update_model_time_step!(model::ModelData)::Nothing
    main_time_loop = model.timestep.parameters.main_time_loop
    timestep_viscoelastic = main_time_loop.timestep_viscoelastic.value
    set_parameter!(model, "timestep", timestep_viscoelastic)
    return nothing
end

function print_loop_info(model_manager::ModelManagerState, ntimestep::Int)::Nothing
    time1 = model_manager.loop_time_start
    time2 = time()
    print_info("Completed time step $ntimestep: cpu(s) $(round(time2-time1, digits=4))", level=1)
    return nothing
end

function print_initial_sealevel(model::ModelData)::Nothing
    sealevel = model.topography.parameters.sealevel
    y_water_ini = sealevel.y_water_ini.value
    y_sealevel = sealevel.y_sealevel.value
    print_info("Initial water level y-coordinate (m) : $y_water_ini", level=2)
    print_info("Current sealevel y-coordinate (m) : $y_sealevel", level=2)
    return nothing
end

""" Time loop for solving transient thermo-mechanical equations.
"""
function run_loop!(model_manager::ModelManagerState)::Nothing
    
    model = model_manager.model
    ntimestep_max = TimeStep.get_maximum_number_of_time_steps(model)
    _ntimestep_initial = TimeStep.get_current_number_of_time_steps(model)
    print_initial_sealevel(model)
    print_info("Initial time step counter: $_ntimestep_initial", level=2)
    SurfaceProcesses.reset_marker_compaction_properties!(model)
    initialize_next_eruption_time!(model)
    inside_flags = GridFuncs.get_marker_inside_flags(model)
    t1 = time()
    for _ntimestep in _ntimestep_initial:ntimestep_max
        model_manager.loop_time_start = time()
        TimeStep.initialize!(model, _ntimestep)
        execute_pre_solver_steps!(model_manager, inside_flags)
        manage_model_output!(model_manager)
        execute_stokes_continuity_solver_steps!(model_manager, inside_flags)
        execute_heat_solver_steps!(model_manager, inside_flags)
        execute_advection_steps!(model_manager, inside_flags)
        advance_model_time!(model)
        inside_flags = execute_post_solver_steps!(model_manager, inside_flags)
        if Stopping.terminate_loop(model)
            break
        end
        TimeStep.update_output_counter!(model)
        print_loop_info(model_manager, _ntimestep)
    end
    t2 = time()
    print_info(">> Total time taken to run time loop: $(t2 - t1) seconds", level=1)
    return nothing
end

function execute_pre_solver_steps!(
    model_manager::ModelManagerState,
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "# Finished executing pre-solver steps" begin
        model = model_manager.model
        output_dir = model_manager.run_paths["output_dir"]

        boolean_options = get_boolean_options(model_manager)
        BoundaryConditions.update_transient_boundary_conditions!(model)
        MarkerFriction.randomize_marker_initial_friction_angle!(model)
        MeltModel.update_marker_melt_model!(model, inside_flags, output_dir)
        RockProperties.update_marker_rock_properties!(model, inside_flags)
        Serpentinization.update_rock_props_for_serpentinization!(model, inside_flags)
        HydrothermalCirculation.update_rock_properties_for_hydrothermal!(model)
        MeltModel.update_rock_properties_for_melt!(model, inside_flags)
        Interpolation.initialize_bilinear_interpolation!(model, inside_flags)
        Interpolation.interpolate_marker_heat_parameters_to_grid!(model, inside_flags)
        Interpolation.interpolate_temperature_back_to_markers_at_startup!(model, inside_flags)
        Rheology.update_marker_rheology!(model, boolean_options, inside_flags)
        Rheology.update_marker_viscoplastic_viscosity!(model, boolean_options, inside_flags)
        MeltModel.update_marker_viscoplastic_viscosity_for_melting!(model, inside_flags)
        MeltModel.update_marker_melt_damage!(model)
        Interpolation.interpolate_marker_stokes_parameters_to_grid!(model, inside_flags)
        Interpolation.interpolate_basic_viscoplastic_viscosity_to_pressure_grid!(model)
    end
    return nothing
end

function execute_stokes_continuity_solver_steps!(
    model_manager::ModelManagerState,
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "# Finished executing Stokes-continuity solver steps" begin
        model = model_manager.model
        loop_input = ModelManager.get_loop_input_struct(model_manager)
        GlobalPlasticityLoop.run_picard_loop!(model, loop_input, inside_flags)
        StokesContinuitySolver.interpolate_grid_viscoplastic_viscosity_to_markers!(
            model, inside_flags)
        TimeStep.adjust_time_step_using_displacement_limit!(model)
        StokesContinuitySolver.execute_post_stokes_solver_steps!(model, inside_flags)
        StokesContinuitySolver.update_marker_stress_for_grid_and_subgrid_changes!(
            model, inside_flags)
    end
    return nothing
end

function execute_heat_solver_steps!(
    model_manager::ModelManagerState,
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "# Finished executing heat solver steps" begin
        HeatSolver.update_marker_temp_for_grid_and_subgrid_changes!(
            model_manager.model, model_manager.config.solver, inside_flags)
    end
    return nothing
end

function execute_advection_steps!(
    model_manager::ModelManagerState,
    inside_flags::Vector{Int8}
)::Nothing
    model = model_manager.model
    @timeit_memit "# Finished executing advection steps" begin
        Advection.advect_markers_and_rotate_stress_tensor!(model, inside_flags)
        Advection.update_marker_strain!(model, inside_flags)
        HeatSolver.apply_sticky_temperature_correction!(model, inside_flags)
        HeatSolver.manage_plume_injection!(model)
    end
    return nothing
end

function execute_post_solver_steps!(
    model_manager::ModelManagerState,
    inside_flags::Vector{Int8}
)::Vector{Int8}
    model = model_manager.model
    boolean_options = get_boolean_options(model_manager)
    output_dir = model_manager.run_paths["output_dir"]

    @timeit_memit "# Finished executing post-solver steps" begin
        MeltModel.Extraction.update_melt_extraction!(model, inside_flags, output_dir)
        SurfaceProcesses.update_topography!(model, inside_flags)
        SurfaceProcesses.Sealevel.update_sealevel!(model)
        update_next_eruption_time!(model)
        Regrid.regrid!(model, boolean_options)
        BoundaryConditions.manage_marker_recycling!(model, inside_flags)
        inside_flags = GridFuncs.get_marker_inside_flags(model)
        SurfaceProcesses.transform_markers_for_surface_processes!(model)
        SurfaceProcesses.apply_salt_deposition_model!(model)
        SurfaceProcesses.reset_marker_compaction_properties!(model)
        Serpentinization.update_marker_serpentinization!(model, inside_flags)
        StokesContinuitySolver.turn_off_gravity!(model)
        TimeStep.manage_time_step_increase_for_variable_step_models!(model)
    end
    return inside_flags
end

function manage_model_output!(
    model_manager::ModelManagerState
)::Nothing
    @timeit_memit "Finished managing model output" begin
        ModelManager.manage_model_output(model_manager)
    end
    # ModelManager.ModelPlots2dManager.manage_loop_plots!(model_manager.model) # TODO: Add this
    return nothing
end

""" Create boolean options for time loop.

Returns
-------
Dict{String, Bool}
    - Dictionary of boolean options for time loop. The keys are the
    option names and the values are the boolean values.
"""
function get_boolean_options(model_manager::ModelManagerState)::Dict{String, Bool}
    bc_bools = BoundaryConditions.get_boolean_options(model_manager.model)
    grid_bools = StaggeredGrid.get_boolean_options(model_manager.model)
    return Dict{String, Bool}(
        "extending_side_boundaries" => get(bc_bools, :extending_side_boundaries, false),
        "no_yielding_in_mobile_wall" => get(bc_bools, :no_yielding_in_mobile_wall, false),
        "no_yielding_in_plate_extension" => get(bc_bools, :no_yielding_in_plate_extension, false),
        "use_boundary_friction_plasticity_model_sandbox" => get(bc_bools, :use_boundary_friction_plasticity_model_sandbox, false),
        "uniform_grid" => get(grid_bools, :uniform_grid, false),
        "use_marker_plasticity" => GlobalPlasticityLoop.is_global_marker_plasticity_loop(model_manager.model)
    )
end

""" Update total model time in seconds (timesum).

Total model time referred to as timesum is advanced using time step updated
based on maximum marker displacement.

Updated Model Parameter
=======================
model.timestep.parameters.main_time_loop
----------------------------------------
timesum.value: float
    - Total model time in seconds.
"""
function advance_model_time!(model::ModelData)::Nothing
    timestep = model.timestep.parameters.main_time_loop.timestep
    timesum = model.timestep.parameters.main_time_loop.timesum
    timesum.value = timesum.value + timestep.value
    return nothing
end

""" Initialize time of next eruption.

Updated Model Parameters
========================
model.melting.parameters.extrusion
--------------------------------
time_of_next_eruption.value: float
    - Time of next eruption in seconds.
"""
function initialize_next_eruption_time!(model::ModelData)::Nothing
    time_of_next_eruption_myr = model.melting.parameters.extrusion.time_of_next_eruption_myr
    eruption_interval_yr = model.melting.parameters.extrusion.eruption_interval_yr.value
    eruption_interval_myr = eruption_interval_yr/1e6
    time_of_next_eruption_myr.value = eruption_interval_myr
    if model.melting.parameters.extrusion.iuse_eruption_interval.value == 1
        print_info("Next eruption time is $(time_of_next_eruption_myr.value) Myr.", level=2)
    end
    return nothing
end

""" Update time of next eruption.

Updated Model Parameter
-----------------------
model.melting.parameters.extrusion.time_of_next_eruption.value: float
    - Time of next eruption in seconds.
"""
function update_next_eruption_time!(model::ModelData)::Nothing
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    sec_per_myr = model.conversion.parameters.sec_per_Myr.value
    timesum_myr = timesum/sec_per_myr

    time_of_next_eruption_myr = model.melting.parameters.extrusion.time_of_next_eruption_myr
    eruption_interval_yr = model.melting.parameters.extrusion.eruption_interval_yr.value
    eruption_interval_myr = eruption_interval_yr/1e6

    if timesum_myr >= time_of_next_eruption_myr.value
        time_of_next_eruption_myr.value += eruption_interval_myr
        if model.melting.parameters.extrusion.iuse_eruption_interval.value == 1
            print_info("Next eruption time is $(time_of_next_eruption_myr.value) Myr.", level=2)
        end
    end
    return nothing
end

end # module 