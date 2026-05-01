"""
    Model.jl

Module containing the customized model setup and execution functions for the 
hot box melting model.

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

const ROOT_PATH_OUTPUT = "/mnt/extradrive1/earthbox_output/hot_box_melting"

# Get the input parameters object so names can be accessed programmatically
const PARAMS = get_eb_parameters()

# Constants are defined for parameters used in multiple initialization functions
const xsize          = 300_000.0 # meters
const ysize          = 150_000.0 # meters
const marker_spacing = 100.0 # meters
const thick_air      = 10_000.0 # meters

function run_case(;case_name::String="case0")::Nothing
    print_info("Running model for case: $case_name")
    model_output_path = get_model_output_path(case_name, ROOT_PATH_OUTPUT)
    eb = setup_model(model_output_path)
    PRINT_SETTINGS.print_performance = true
    run_time_steps(
        eb,
        ntimestep_max         = 10,
        timestep_viscoelastic = ConversionFuncs.years_to_seconds(50_000.0),
        timestep_out          = ConversionFuncs.years_to_seconds(50_000.0)
    )
    return nothing
end

function setup_model(output_dir::String)::EarthBoxState
    eb = EarthBoxState(
        restart_from_backup         = false,
        xnum                        = 151,
        ynum                        = 76,
        xsize                       = xsize,
        ysize                       = ysize,
        dx_marker                   = marker_spacing,
        dy_marker                   = marker_spacing,
        use_mumps                   = false,
        nprocs                      = 8,
        paths                       = Dict("output_dir" => output_dir)
    )
    initialize_model_input(eb)
    return eb
end

function initialize_model_input(eb::EarthBoxState)::Nothing
    initialize_marker_output(eb)
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
    initialize_melt_model(eb)
    initialize_surface_processes_model(eb)
    return nothing
end

function initialize_marker_output(eb::EarthBoxState)::Nothing
    keys = eb.model_manager.config.output.marker_output_keys
    marker_output = eb.model_manager.config.output.marker_output
    
    marker_output[keys.marker_rho]                  = false
    marker_output[keys.marker_serpentinization]     = false
    marker_output[keys.marker_extractable_meltfrac] = false
    marker_output[keys.marker_extracted_meltfrac]   = false
    marker_output[keys.marker_TK]                   = true
    
    return nothing
end

function initialize_staggered_grid(eb::EarthBoxState)::Nothing
    StaggeredGrid.initialize!(
        eb.model_manager.model, grid_type = StaggeredGrid.option_names.UniformGrid)
    return nothing
end

function initialize_geometry(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    MaterialGeometry.StickyAirGeometry.initialize!(model, thick_air = thick_air)
    return nothing
end

function initialize_boundary_conditions(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    BoundaryConditions.initialize!(
        model, model_type = bc_model_type_names.BoxConvection)
    BoundaryConditions.Pressure.initialize!(model, pressure_bc=1e5)
    BoundaryConditions.Temperature.initialize!(
        model,
        temperature_top = ConversionFuncs.celsius_to_kelvin(0.0),
        temperature_bottom = ConversionFuncs.celsius_to_kelvin(1330.0)
        )
    return nothing
end
 
function initialize_marker_coordinates(eb::EarthBoxState)::Nothing
    Markers.MarkerCoordinates.initialize!(
        eb.model_manager.model, marker_distribution=marker_distribution_names.Regular)
    return nothing
end

function initialize_marker_materials(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    Markers.MarkerMaterials.initialize!(
        model,
        material_model                = material_model_names.SimpleSediment,
        paths                         = Dict("materials_library_file" => MATERIAL_COLLECTION.path),
        materials_input_dict          = get_materials_input_dict(),
        viscosity_min                 = 1e18, # Pa.s
        viscosity_max                 = 1e26, # Pa.s
    )
    return nothing
end

function initialize_marker_temperature(eb::EarthBoxState)::Nothing
    Markers.MarkerTemperature.initialize!(
        eb.model_manager.model,
        initial_temperature_model=initial_temperature_names.HotBox
    )
    return nothing
end

function initialize_marker_plasticity(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    Markers.MarkerFriction.initialize!(model, initialization_model = 0)
    Markers.MarkerCohesion.initialize!(model)
    return nothing
end

function initialize_rock_properties(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    RockProperties.initialize!(
        model,
        density_model              = density_model_names.Liao14,
        thermal_conductivity_model = thermal_conductivity_model_names.Liao14,
        rhocp_model                = rhocp_model_names.Constant
    )
    return nothing
end

function initialize_stokes_continuity_solver(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    StokesContinuitySolver.initialize!(
        model,
        velocity_type                = velocity_type_names.VelocityFromStokesSolver,
        gravity_y                    = 9.8, # m/s^2
        iuse_interface_stabilization = 1,
    )
    GlobalPlasticityLoop.initialize!(
        model,
        global_plasticity_loop = global_plasticity_names.NodalPlasticityLoop,
        tolerance_picard       = 0.01,
        nglobal                = 2,
    )
    return nothing
end

function initialize_heat_solver(eb::EarthBoxState)::Nothing
    HeatSolver.initialize!(
        eb.model_manager.model,
        iuse_heat              = 1,
        iuse_adiabatic_heating = 1,
        iuse_shear_heating     = 1,
        max_temp_change        = 70.0, # K
    )
    return nothing
end

function initialize_advection_model(eb::EarthBoxState)::Nothing
    Advection.initialize!(
        eb.model_manager.model,
        advection_scheme         = advection_scheme_names.RungeKutta4thOrder,
        marker_cell_displ_max    = 1.0, # fraction
        subgrid_diff_coef_temp   = 1.0,
        subgrid_diff_coef_stress = 1.0,
    )    
    return nothing
end

function initialize_interpolation_model(eb::EarthBoxState)::Nothing
    Interpolation.initialize!(
        eb.model_manager.model,
        iuse_harmonic_avg_normal_viscosity = 1,
        iuse_initial_temp_interp           = 1,
    )
    return nothing
end

function initialize_melt_model(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    MeltModel.initialize!(
        model,
        iuse_melting                         = 1,
        iuse_melt_viscosity                  = 1,
        iuse_melt_thermal_props              = 1,
        viscosity_melt                       = 1e18, # Pa.s
        iuse_depletion_density               = 0,
        iuse_exponential_viscosity_reduction = 0,
    )
    MeltModel.Extraction.initialize!(
        model,
        iuse_extraction                 = 1,
        extraction_fraction             = 0.99, # fraction
        iuse_shallow_mantle_injection   = 1,
        iuse_random_injection_subdomain = 1,
        iuse_normal_injection_subdomain = 1,
        smoothing_radius_drainage       = 10_000.0, # m
        smoothing_radius_fractionation  = 10_000.0, # m
        characteristic_injection_width  = 10_000.0, # m
        mantle_search_width             = 300_000.0, # m
        number_of_injection_subdomains  = 10,
        magma_height_limit              = 3000.0, # m
        emplacement_temperature         = ConversionFuncs.celsius_to_kelvin(1_200.0),
        iuse_gabbroic_fractionation     = 1,
        fractionation_threshold_limit   = 2000.0, # m
        maximum_shallow_injection_depth = 40_000.0, # m
    )
    MeltModel.Extrusion.initialize!(
        model,
        iuse_extrusion                           = 1,
        extrusion_volume_factor                  = 0.2, # fraction
        extrusion_volume_factor_max              = 0.2, # fraction
        characteristic_magmatic_crust_height_min = 6_000.0, # m
        characteristic_magmatic_crust_height_max = 6_000.0, # m
        width_eruption_domain_fixed              = 2_500.0, # m
        width_eruption_domain_fixed_max          = 2_500.0, # m
        characteristic_flow_length_subaerial     = 100_000.0, # m
        characteristic_flow_length_submarine     = 10_000.0, # m
        residual_lava_thickness_subaerial        = 30.0, # m
        residual_lava_thickness_submarine        = 60.0, # m
        iuse_random_eruption_location            = 1,
        iuse_normal_eruption_location            = 1,
        porosity_initial_lava_flow               = 0.0, # fraction
        decay_depth_lava_flow                    = 2000.0, # m
        decimation_factor                        = 4, # None
        iuse_eruption_interval                   = 0, # None
        eruption_interval_yr                     = 50_000.0, # yr
    )
    return nothing
end

function initialize_surface_processes_model(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    SurfaceProcesses.Topography.initialize!(
        model,
        iuse_topo            = 1,
        node_advection       = topo_node_advection_names.RungeKuttaWithInterp,
        dx_topo              = 200.0, # m
        topo_xsize           = xsize, # m
        nsmooth_top_bottom   = 1,
        marker_search_factor = 2.0,
    )
    SurfaceProcesses.Sealevel.initialize!(
        model,
        option_name               = sealevel_option_names.Constant,
        y_sealevel                = thick_air, # m
    )
    SurfaceProcesses.SedimentTransport.initialize!(
        model,
        iuse_downhill_diffusion         = 1,
        subaerial_transport_coefficient = 1.0e-2,
        subaerial_slope_diffusivity     = ConversionFuncs.m2_yr_to_m2_s(0.25),
        submarine_slope_diffusivity     = ConversionFuncs.m2_yr_to_m2_s(100.0),
        submarine_diffusion_decay_depth = 2000.0, # m
        precipitation_rate              = ConversionFuncs.m_yr_to_m_s(1.0),
        pelagic_sedimentation_rate      = 0.0, # mm/yr
        iuse_compaction_correction      = 1,
        transport_timestep              = ConversionFuncs.years_to_seconds(10_000.0)
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
