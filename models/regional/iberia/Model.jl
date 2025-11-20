"""
    Model.jl

Module containing the customized model setup and execution functions for the 
regional/iberia model.

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

const ROOT_PATH_OUTPUT = "/mnt/extradrive1/earthbox_output/regional/iberia"

# Get the input parameters object so names can be accessed programmatically
const PARAMS = get_eb_parameters()

# Constants are defined for parameters used in multiple initialization functions
const adiabatic_gradient       = 0.4 # K/km
const xsize                    = 1_000_000.0 # meters
const ysize                    = 160_000.0 # meters
const thick_air                = 10_000.0 # meters
const thick_crust              = 35_000.0 # meters
const thick_upper_crust        = 22_000.0 # meters
const thick_lith               = 125_000.0 # meters
const marker_spacing           = 100.0 # meters (50 m for high resolution case)
const grid_spacing_high_res    = 500.0 # meters (200 m for high resolution case)
const avg_grid_spacing_low_res = 2000.0 # meters
const temperature_base_lith_celsius = 1345.0 # C, this is the final cooler temperature at the base of the lithosphere

function run_case(;case_name::String="case0")::Nothing
    print_info("Running model for case: $case_name")
    model_output_path = get_model_output_path(case_name, ROOT_PATH_OUTPUT)
    eb = setup_model(model_output_path)
    PRINT_SETTINGS.print_performance = true
    run_time_steps(
        eb,
        make_backup           = true,
        ntimestep_max         = 2500,
        timestep_viscoelastic = ConversionFuncs.years_to_seconds(50_000.0),
        timestep_out          = ConversionFuncs.years_to_seconds(200_000.0)
    )
    return nothing
end

function setup_model(output_dir::String)::EarthBoxState
    ttype_refinement_parameters = Dict(
        PARAMS.xo_highres.name => 350_000.0,
        PARAMS.xf_highres.name => 650_000.0,
        PARAMS.yf_highres.name => 70_000.0,
        PARAMS.dx_highres.name => grid_spacing_high_res,
        PARAMS.dx_lowres.name  => avg_grid_spacing_low_res,
        PARAMS.dy_highres.name => grid_spacing_high_res,
        PARAMS.dy_lowres.name  => avg_grid_spacing_low_res
    )
    eb = EarthBoxState(
        restart_from_backup         = false,
        xsize                       = xsize,
        ysize                       = ysize,
        dx_marker                   = marker_spacing,
        dy_marker                   = marker_spacing,
        ttype_refinement_parameters = ttype_refinement_parameters,
        use_mumps                   = true,
        nprocs                      = 8,
        analysis_method             = :PARALLEL,
        parallel_ordering_method    = :PTSCOTCH,
        paths                       = Dict("output_dir" => output_dir)
    )
    initialize_model_input(eb)
    return eb
end

function initialize_model_input(
    eb::EarthBoxState
)::Nothing
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
    
    marker_output[keys.marker_rho]                  = true
    marker_output[keys.marker_serpentinization]     = true
    marker_output[keys.marker_extractable_meltfrac] = false
    marker_output[keys.marker_extracted_meltfrac]   = false
    
    return nothing
end

function initialize_staggered_grid(eb::EarthBoxState)::Nothing
    StaggeredGrid.initialize!(
        eb.model_manager.model, grid_type = StaggeredGrid.option_names.TtypeRefinedGrid)
    return nothing
end

function initialize_geometry(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    MaterialGeometry.StickyAirGeometry.initialize!(model, thick_air = thick_air)
    MaterialGeometry.EarthLayering.initialize!(
        model,
        thick_lith        = thick_lith,
        thick_crust       = thick_crust,
        thick_upper_crust = thick_upper_crust,
    )
    MaterialGeometry.LithoStrongZones.initialize!(model, iuse_strong_zones = 1)
    return nothing
end

function initialize_boundary_conditions(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    BoundaryConditions.initialize!(
        model, model_type = bc_model_type_names.LithosphericExtensionFixedBoundaries)
    BoundaryConditions.Pressure.initialize!(model, pressure_bc=1e5)
    BoundaryConditions.Temperature.initialize!(
        model, 
        temperature_top = ConversionFuncs.celsius_to_kelvin(0.0),
        temperature_bottom = ConversionFuncs.celsius_to_kelvin(
              temperature_base_lith_celsius 
            + adiabatic_gradient*(ysize - thick_air - thick_lith)/1000.0
            )
        )
    BoundaryConditions.Velocity.initialize!(
        model, full_velocity_extension = ConversionFuncs.cm_yr_to_m_s(0.14))
    BoundaryConditions.VelocityStep.initialize!(
        model,
        iuse_velocity_step                = 1,
        velocity_step_factor              = 0.4/0.14, # 0.14 cm/yr to 0.4 cm/yr
        timestep_adjustment_factor        = 40_000.0/50_000.0, # 50_000.0 to 40_000.0 yr
        velocity_step_time                = 50.0, # Myr
        velocity_second_step_factor       = 1.0/0.4, # 0.4 cm/yr to 1 cm/yr
        timestep_second_adjustment_factor = 15_000.0/40_000.0, # 40_000 yr to 15_000 yr
        velocity_second_step_time         = 90.0
    )
    return nothing
end
 
function initialize_marker_coordinates(eb::EarthBoxState)::Nothing
    Markers.MarkerCoordinates.initialize!(
        eb.model_manager.model, marker_distribution=marker_distribution_names.Randomized)
    return nothing
end

function initialize_marker_materials(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    Markers.MarkerMaterials.initialize!(
        model,
        material_model                = material_model_names.LithosphericExtensionLateralStrongZones,
        paths                         = Dict("materials_library_file" => MATERIAL_COLLECTION.path),
        materials_input_dict          = get_materials_input_dict(),
        viscosity_min                 = 1e18, # Pa.s
        viscosity_max                 = 1e26, # Pa.s
        iuse_fluid_pressure_for_yield = 1,
        plastic_healing_rate          = 0.0, # 1/s
    )
    Markers.MarkerMaterials.MarkerStressLimits.initialize!(
        model,
        powerlaw_stress_min = 1e4, # Pa
        yield_stress_min    = 1e6, # Pa
        yield_stress_max    = 1e32, # Pa
    )
    Markers.MarkerMaterials.MarkerViscousStrainSoftening.initialize!(
        model, 
        iuse_viscous_strain_soft = 0, 
        vsoftfac                 = 30.0,
        )
    return nothing
end

function initialize_marker_temperature(eb::EarthBoxState)::Nothing
    Markers.MarkerTemperature.initialize!(
        eb.model_manager.model,
        initial_temperature_model=initial_temperature_names.AnalyticalThreeLayer,
        parameters=Dict(
            PARAMS.amplitude_perturbation.name      => 5_000.0,
            PARAMS.width_perturbation.name          => 10_000.0,
            PARAMS.temperature_base_lith.name       => ConversionFuncs.celsius_to_kelvin(temperature_base_lith_celsius),
            PARAMS.thick_thermal_lithosphere.name   => thick_lith,
            PARAMS.adiabatic_gradient.name          => adiabatic_gradient,
            PARAMS.conductivity_upper_crust.name    => 2.25, # W/m/K
            PARAMS.conductivity_lower_crust.name    => 2.0, # W/m/K
            PARAMS.conductivity_mantle.name         => 2.0, # W/m/K
            PARAMS.heat_production_upper_crust.name => 1.8e-6, # W/m^3
            PARAMS.heat_production_lower_crust.name => 0.5e-6, # W/m^3
            PARAMS.heat_production_mantle.name      => 0.0, # W/m^3
        )
    )
    return nothing
end

function initialize_marker_plasticity(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    Markers.MarkerFriction.initialize!(
        model,
        initialization_model  = 0, # regular initialization of friction coefficients
        iuse_random_fric_time = 1, # Activate time-dependent randomization of friction  coefficients
        randomization_factor  = 10.0, # Time-dependent randomization factor for friction coefficients
    )
    Markers.MarkerCohesion.initialize!(model)
    Markers.MarkerPreexponential.initialize!(model)
    return nothing
end

function initialize_rock_properties(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    RockProperties.initialize!(
        model,
        density_model              = density_model_names.Liao14,
        thermal_conductivity_model = thermal_conductivity_model_names.SekiguchiWaples,
        rhocp_model                = rhocp_model_names.TemperatureDependentWaples,
        iuse_melt_lens             = 0,
        iuse_sed_porosity          = 1,
    )
    HydrothermalCirculation.initialize!(
        model,
        iuse_hydrothermal                          = 1,
        hydrothermal_smoothing_factor              = 0.75,
        hydrothermal_max_depth                     = 4000.0, # m
        hydrothermal_max_temperature               = 600.0, # C
        hydrothermal_nusselt_number                = 2.0,
        iuse_plastic_strain_rate_for_hydrothermal  = 1,
        hydrothermal_plastic_strain_rate_reference = 1e-14, # 1/s
        iuse_plastic_strain_for_hydrothermal       = 1,
        hydrothermal_plastic_strain_reference      = 0.5,
        hydrothermal_decay_length                  = 25_000.0, # m
        hydrothermal_buffer_distance               = 25_000.0, # m
        sediment_thickness_threshold               = 2500.0, # m
    )
    Serpentinization.initialize!(
        model,
        iuse_serpentinization                = 1, # Yes
        serpentinization_temperature         = 340.5, # C
        maximum_serpentinization_depth       = 4_000.0, # m
        maximum_serpentinization_rate        = 1e-11, # 1/s
        nominal_strain_rate_serpentinization = 1e-13, # 1/s
    )
    return nothing
end

function initialize_stokes_continuity_solver(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    StokesContinuitySolver.initialize!(
        model,
        velocity_type                = velocity_type_names.VelocityFromStokesSolver,
        gravity_y                    = 9.8, # m/s^2
        iuse_interface_stabilization = 1
    )
    GlobalPlasticityLoop.initialize!(
        model,
        global_plasticity_loop = global_plasticity_names.NodalPlasticityLoop,
        tolerance_picard       = 1e-2,
        nglobal                = 3,
    )
    return nothing
end

function initialize_heat_solver(eb::EarthBoxState)::Nothing
    HeatSolver.initialize!(
        eb.model_manager.model,
        iuse_heat              = 1,
        iuse_adiabatic_heating = 1,
        iuse_shear_heating     = 1,
        iuse_sticky_correction = 1,
        max_temp_change        = 70.0, # K
    )
    return nothing
end

function initialize_advection_model(eb::EarthBoxState)::Nothing
    Advection.initialize!(
        eb.model_manager.model,
        advection_scheme                  = advection_scheme_names.RungeKutta4thOrder,
        #iuse_local_adaptive_time_stepping = 1,
        marker_cell_displ_max             = 1.0, # fraction
        subgrid_diff_coef_temp            = 1.0,
        subgrid_diff_coef_stress          = 1.0,
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
        iuse_depletion_density               = 1,
        iuse_exponential_viscosity_reduction = 1,
    )
    MeltModel.MeltDamage.initialize!(
        model,
        iuse_melt_damage                   = 1,
        melt_damage_distance               = 2_500.0, # m
        melt_damage_factor                 = 10.0, # None
        iuse_probabilistic_melt_damage     = 1,
        maximum_damage_probability         = 0.8, # fraction
        intermediate_damage_probability    = 0.1, # fraction
        magmatic_crust_height_threshold    = 500.0, # m
        magmatic_crust_height_minimum      = 750.0, # m
        magmatic_crust_height_intermediate = 2_000.0, # m
        magmatic_crust_height_maximum      = 3_000.0, # m
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
        mantle_search_width             = 200_000.0, # m
        number_of_injection_subdomains  = 10, # None
        magma_height_limit              = 1000.0, # m
        emplacement_temperature         = ConversionFuncs.celsius_to_kelvin(1_200.0),
        iuse_gabbroic_fractionation     = 1,
        fractionation_threshold_limit   = 2000.0, # m
        maximum_shallow_injection_depth = 25_000.0, # m
    )
    MeltModel.Extrusion.initialize!(
        model,
        iuse_extrusion                           = 1,
        extrusion_volume_factor                  = 0.06, # fraction
        extrusion_volume_factor_max              = 0.5, # fraction
        characteristic_magmatic_crust_height_min = 6_000.0, # m
        characteristic_magmatic_crust_height_max = 7_500.0, # m
        width_eruption_domain_fixed              = 2_500.0, # m
        width_eruption_domain_fixed_max          = 2_500.0, # m
        characteristic_flow_length_subaerial     = 20_000.0, # m
        characteristic_flow_length_submarine     = 2_000.0, # m
        residual_lava_thickness_subaerial        = 30.0, # m
        residual_lava_thickness_submarine        = 30.0, # m
        iuse_random_eruption_location            = 1,
        iuse_normal_eruption_location            = 1,
        porosity_initial_lava_flow               = 0.0, # fraction
        decay_depth_lava_flow                    = 2000.0, # m
        decimation_factor                        = 4,
        iuse_eruption_interval                   = 1,
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
        nsmooth_top_bottom   = 2,
        marker_search_factor = 2.0,
    )
    SurfaceProcesses.Sealevel.initialize!(
        model,
        option_name               = sealevel_option_names.AveragePressure,
        y_sealevel                = thick_air, # m
        base_level_shift_end_time = 0.0, # Myr
        base_level_shift          = 0.0, # m
    )
    SurfaceProcesses.Sealevel.RelativeBaseLevel.initialize!(
        model,
        iuse_linear_segments                  = 1,
        thickness_upper_continental_crust_ref = thick_upper_crust, # m
        thickness_lower_continental_crust_ref = thick_crust - thick_upper_crust, # m
        thickness_lithosphere_ref             = thick_lith, # m
        gridy_spacing_ref                     = 100.0, # m
        temperature_top_ref                   = 0.0, # C
        temperature_moho_ref                  = 600.0, # C
        temperature_base_lith_ref             = temperature_base_lith_celsius, # C
        adiabatic_gradient_ref                = adiabatic_gradient, # K/km
    )
    SurfaceProcesses.SedimentTransport.initialize!(
        model,
        iuse_downhill_diffusion         = 1,
        subaerial_transport_coefficient = 1.0e-4, # None
        subaerial_slope_diffusivity     = ConversionFuncs.m2_yr_to_m2_s(0.25),
        submarine_slope_diffusivity     = ConversionFuncs.m2_yr_to_m2_s(100.0),
        submarine_diffusion_decay_depth = 1000.0, # m
        precipitation_rate              = ConversionFuncs.m_yr_to_m_s(1.0),
        pelagic_sedimentation_rate      = 0.0, # mm/yr
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
