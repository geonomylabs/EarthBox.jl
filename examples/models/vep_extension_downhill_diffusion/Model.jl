"""
    Model.jl

Module containing the customized model setup and execution functions for the 
visco-elasto-plastic extension with downhill diffusion model.

Model output is sent to the model output directory named `<case_name>_output` 
where `<case_name>` is the name of the model case (e.g. "case1"). The model 
output path is defined using the `case_name` and the `ROOT_PATH_OUTPUT`constant 
unless the `model_output_path` is provided via command line argument.

See [Building and Executing Models with the EarthBox API](@ref) in the EarthBox
docs for more details on command line arguments and usage.

Quick Start:

From the REPL, run:

```julia
include("Model.jl")
Model.run_case(case_name="case1")
```

If no `<case_name>` is provided, then the default case name is `case0`.

From command-line, run:

```bash
julia Model.jl case_name=<case_name>
```
"""
module Model

using EarthBox
include("Materials.jl")
import .Materials: get_materials_input_dict, MATERIAL_COLLECTION

const ROOT_PATH_OUTPUT = "/mnt/extradrive1/earthbox_output/vep_extension_downhill_diffusion"

# Get the input parameters object so names can be accessed programmatically
const PARAMS = get_eb_parameters()

# Constants are defined for parameters used in multiple initialization functions
const xsize                    = 500_000.0 # meters
const ysize                    = 140_000.0 # meters
const thick_air                = 10_000.0 # meters
const thick_crust              = 32_000.0 # meters
const thick_upper_crust        = 22_000.0 # meters
const thick_lith               = 122_000.0 # meters
const marker_spacing           = 100.0 # meters (50 m for high resolution case)
const grid_spacing_high_res    = 500.0 # meters (200 m for high resolution case)
const avg_grid_spacing_low_res = 2000.0 # meters

function run_case(;case_name::String="case0")::Nothing
    print_info("Running model for case: $case_name")
    model_output_path = get_model_output_path(case_name, ROOT_PATH_OUTPUT)
    eb = setup_model(model_output_path=model_output_path)
    PRINT_SETTINGS.print_performance = true
    run_time_steps(
        eb,
        make_backup           = true,
        ntimestep_max         = 600,
        timestep_viscoelastic = ConversionFuncs.years_to_seconds(50_000.0),
        timestep_out          = ConversionFuncs.years_to_seconds(500_000.0)
    )
    return nothing
end

function setup_model(;model_output_path::String)::EarthBoxState
    ttype_refinement_parameters = Dict(
        PARAMS.xo_highres.name => 150_000.0, # meters
        PARAMS.xf_highres.name => 350_000.0, # meters
        PARAMS.yf_highres.name => 50_000.0, # meters
        PARAMS.dx_highres.name => grid_spacing_high_res, # meters (200 m for high resolution case)
        PARAMS.dx_lowres.name  => avg_grid_spacing_low_res, # meters
        PARAMS.dy_highres.name => grid_spacing_high_res, # meters (200 m for high resolution case)
        PARAMS.dy_lowres.name  => avg_grid_spacing_low_res # meters
    )
    eb = EarthBoxState(
        restart_from_backup         = false,
        xsize                       = xsize, # meters
        ysize                       = ysize, # meters
        dx_marker                   = marker_spacing, # meters (50 m for high resolution case)
        dy_marker                   = marker_spacing, # meters (50 m for high resolution case)
        ttype_refinement_parameters = ttype_refinement_parameters,
        use_mumps                   = true,
        nprocs                      = 8,
        paths                       = Dict("output_dir" => model_output_path)
    )
    initialize_model_input(eb)
    return eb
end

function initialize_model_input(eb::EarthBoxState)::Nothing
    initialize_staggered_grid(eb)
    initialize_geometry(eb)
    initialize_boundary_conditions(eb)
    initialize_marker_coordinates(eb)
    initialize_marker_materials(eb)
    initialize_marker_temperature(eb)
    initialize_marker_plasticity(eb)
    initialize_rock_properties(eb)
    initialize_stokes_continuity_solver(eb)
    initialize_heat_solver(eb)
    initialize_advection_model(eb)
    initialize_interpolation_model(eb)
    initialize_surface_processes_model(eb)
    return nothing
end

function initialize_staggered_grid(eb::EarthBoxState)::Nothing
    StaggeredGrid.initialize!(
        eb.model_manager.model, grid_type = grid_type_names.TtypeRefinedGrid)
    return nothing
end

function initialize_geometry(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    MaterialGeometry.StickyAirGeometry.initialize!(model, thick_air=thick_air)
    MaterialGeometry.EarthLayering.initialize!(
        model,
        thick_lith        = thick_lith, # meters
        thick_crust       = thick_crust, # meters
        thick_upper_crust = thick_upper_crust, # meters
    )
    MaterialGeometry.WeakFault.initialize!(
        model,
        fault_dip_degrees = 45.0, # degrees
        fault_thickness = 400.0, # meters
        x_initial_fault = 250_000.0, # meters
        fault_height = 40_000.0, # meters
        )
    return nothing
end

function initialize_boundary_conditions(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    BoundaryConditions.initialize!(
        model, model_type = bc_model_type_names.LithosphericExtensionFixedBoundaries)
    BoundaryConditions.Pressure.initialize!(model, pressure_bc=1e5)
    BoundaryConditions.Velocity.initialize!(
        model, full_velocity_extension = ConversionFuncs.cm_yr_to_m_s(0.1))
    BoundaryConditions.VelocityStep.initialize!(
        model,
        iuse_velocity_step = 1,
        velocity_step_factor = 10.0,
        timestep_adjustment_factor = 0.5,
        velocity_step_time = 0.5, # Myr
    )
    return nothing
end
 
function initialize_marker_coordinates(eb::EarthBoxState)::Nothing
    Markers.MarkerCoordinates.initialize!(
        eb.model_manager.model, 
        marker_distribution=marker_distribution_names.Randomized
        )
    return nothing
end

function initialize_marker_materials(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    Markers.MarkerMaterials.initialize!(
        model,
        material_model       = material_model_names.LithosphericExtensionWeakFault,
        paths                = Dict("materials_library_file" => MATERIAL_COLLECTION.path),
        materials_input_dict = get_materials_input_dict(),
        viscosity_min        = 1e18, # Pa.s
        viscosity_max        = 1e26 # Pa.s
    )
    Markers.MarkerMaterials.MarkerStressLimits.initialize!(
        model, 
        yield_stress_min    = 0.0, # Pa
        yield_stress_max    = 1.0e9 # Pa
    )
    return nothing
end

function initialize_marker_temperature(eb::EarthBoxState)::Nothing
    Markers.MarkerTemperature.initialize!(
        eb.model_manager.model,
        initial_temperature_model=initial_temperature_names.Uniform,
        temperature_uniform_kelvins=ConversionFuncs.celsius_to_kelvin(0.0)
    )
    return nothing
end

function initialize_marker_plasticity(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    Markers.MarkerFriction.initialize!(
        model, initialization_model = friction_init_names.Regular)
    Markers.MarkerCohesion.initialize!(model)
    Markers.MarkerPreexponential.initialize!(model)
    return nothing
end

function initialize_rock_properties(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    RockProperties.initialize!(
        model,
        density_model              = density_model_names.Liao14,
        thermal_conductivity_model = thermal_conductivity_model_names.Liao14,
        iuse_sed_porosity          = 1
    )
    return nothing
end

function initialize_stokes_continuity_solver(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    StokesContinuitySolver.initialize!(
        model,
        velocity_type                = velocity_type_names.VelocityFromStokesSolver,
        gravity_y                    = 9.8, # m/s^2
        iuse_interface_stabilization = 1 # avoid the drunken sailor problem
    )
    GlobalPlasticityLoop.initialize!(
        model,
        global_plasticity_loop = global_plasticity_names.NodalPlasticityLoop,
        tolerance_picard       = 1e-2,
        nglobal                = 3
    )
    return nothing
end

function initialize_heat_solver(eb::EarthBoxState)::Nothing
    HeatSolver.initialize!(eb.model_manager.model, iuse_heat = 0)
    return nothing
end

function initialize_advection_model(eb::EarthBoxState)::Nothing
    Advection.initialize!(
        eb.model_manager.model,
        advection_scheme                  = advection_scheme_names.RungeKutta4thOrder,
        marker_cell_displ_max             = 0.5, # fraction
        subgrid_diff_coef_stress          = 1.0,
    )    
    return nothing
end

function initialize_interpolation_model(eb::EarthBoxState)::Nothing
    Interpolation.initialize!(
        eb.model_manager.model, iuse_harmonic_avg_normal_viscosity = 1)
    return nothing
end

function initialize_surface_processes_model(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    SurfaceProcesses.Topography.initialize!(
        model,
        iuse_topo            = 1,
        node_advection       = topo_node_advection_names.RungeKuttaWithInterp,
        dx_topo              = 200.0, # meters
        topo_xsize           = xsize, # meters
        nsmooth_top_bottom   = 2,
        marker_search_factor = 4.0    
    )
    SurfaceProcesses.Sealevel.initialize!(
        model,
        option_name               = sealevel_option_names.AveragePressure,
        y_sealevel                = thick_air, # meters
        base_level_shift_end_time = 1000.0, # Myr
        base_level_shift          = 0.0, # meters
    )
    SurfaceProcesses.Sealevel.RelativeBaseLevel.initialize!(
        model,
        iuse_linear_segments                  = 1,
        thickness_upper_continental_crust_ref = thick_upper_crust, # meters
        thickness_lower_continental_crust_ref = thick_crust - thick_upper_crust, # meters
        thickness_lithosphere_ref             = thick_lith, # meters
        gridy_spacing_ref                     = 100.0, # meters
        temperature_top_ref                   = 0.0, # Celsius
        temperature_moho_ref                  = 600.0, # Celsius
        temperature_base_lith_ref             = 1330.0, # Celsius
        adiabatic_gradient_ref                = 0.4, # K/km
    )
    SurfaceProcesses.SedimentTransport.initialize!(
        model,
        iuse_downhill_diffusion         = 1,
        subaerial_transport_coefficient = 1.0e-2,
        subaerial_slope_diffusivity     = ConversionFuncs.m2_yr_to_m2_s(0.25),
        submarine_slope_diffusivity     = ConversionFuncs.m2_yr_to_m2_s(100.0),
        submarine_diffusion_decay_depth = 1000.0, # m
        precipitation_rate              = ConversionFuncs.m_yr_to_m_s(1.0),
        pelagic_sedimentation_rate      = 0.0, # m/s
        iuse_compaction_correction      = 1
    )
    return nothing
end

function main()::Nothing
    run_case(case_name=get_case_name_from_cl_args())
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Model.main()
end
