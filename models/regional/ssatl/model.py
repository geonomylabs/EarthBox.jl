#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Script for working with models where parameters are input through the API.

Usage
-----
To run model:
    % python model.py model_case_name

"""
from paths import get_output_path
from materials import get_materials_input_dict
from plot import get_plot_parameters
from case_inputs import get_case_inputs
from earthbox import EarthBox
from earthbox.plot_tools.utils.loop_plots import setup_loop_plotting
from earthbox.plot_tools import MarkerPlotNames
from earthbox.utils.get_paths import get_model_case_name
from earthbox.grids import TtypeGrid
from earthbox.input_tools.case_inputs_tools import CaseInputsKeys
from earthbox.input_tools import input_preprocessing


def run_model(
        use_loop_plotting: bool=False,
        model_case_name: str='case0'
) -> None:
    """ Set up model and run model time steps
    """
    output_dir = get_output_path(
        model_case_name=model_case_name, drive_number_id=1)
    storage_dir = get_output_path(
        model_case_name=model_case_name, drive_number_id=2)
    case_inputs = get_case_inputs(model_case_name)

    keys = CaseInputsKeys()
    eb = setup_model(case_inputs, output_dir, storage_dir)

    eb.model_manager.config.solver.use_mumps = case_inputs['use_mumps']
    eb.model_manager.config.solver.nproc = 8

    if use_loop_plotting:
        setup_loop_plotting(
            eb, activate_loop_plots, get_plot_parameters, model_case_name)

    (
        time_step_parameters
    ) = input_preprocessing.get_time_step_parameters_for_velocity_stepping(case_inputs)

    years_to_seconds = eb.units.years_to_seconds
    eb.run_time_steps(
        number_of_time_steps=time_step_parameters['total_number_of_time_steps'],
        viscoelastic_time_step_seconds=years_to_seconds(
            time_step_parameters['viscoelastic_time_step_yrs_initial']),
        output_time_step_seconds=years_to_seconds(
            case_inputs[keys.output_time_step_years]),
        make_backup=True
        )


def setup_model(
        case_inputs: dict[str | float | int],
        output_dir: str,
        storage_dir: str | None = None
) -> EarthBox:
    """ Setup model
    """
    keys = CaseInputsKeys()
    ttype_grid = TtypeGrid(
        model_width_meters=case_inputs[keys.model_width_meters],
        model_height_meters=case_inputs[keys.model_height_meters],
        xo_highres=case_inputs[keys.xo_highres],
        xf_highres=case_inputs[keys.xf_highres],
        dx_highres=case_inputs[keys.dx_highres],
        dx_lowres=case_inputs[keys.dx_lowres],
        yf_highres=case_inputs[keys.yf_highres],
        dy_highres=case_inputs[keys.dy_highres],
        dy_lowres=case_inputs[keys.dy_lowres],
        dx_marker=case_inputs[keys.dx_marker],
        dy_marker=case_inputs[keys.dy_marker],
        dx_topo=case_inputs[keys.dx_topo]
        )
    eb = EarthBox(
        restart_from_backup=False,
        paths={'output_dir': output_dir, 'storage_dir': storage_dir},
        xnum=ttype_grid.calculate_number_of_basic_nodes_x(),
        ynum=ttype_grid.calculate_number_of_basic_nodes_y(),
        xsize=ttype_grid.model_width_meters,
        ysize=ttype_grid.model_height_meters,
        nmarkers_cell_x=ttype_grid.calculate_number_of_markers_per_cell_x(),
        nmarkers_cell_y=ttype_grid.calculate_number_of_markers_per_cell_y(),
        )
    initialize_model_input(eb, case_inputs, ttype_grid)
    return eb


def initialize_model_input(
        eb: EarthBox,
        case_inputs: dict[str | float | int],
        ttype_grid: TtypeGrid
) -> None:
    """ User defined functions to initialize input via the EarthBox API.
    """
    initialize_marker_output(eb)
    initialize_staggered_grid(eb, ttype_grid)
    initialize_boundary_conditions(eb, case_inputs)
    initialize_marker_coordinates(eb)
    initialize_geometry(eb, case_inputs)
    initialize_marker_materials(eb, case_inputs)
    initialize_marker_temperature(eb, case_inputs)
    initialize_marker_plasticity(eb)
    initialize_rock_properties(eb, case_inputs)
    initialize_stokes_continuity_solver(eb, case_inputs)
    initialize_heat_solver(eb)
    initialize_advection_model(eb, case_inputs)
    initialize_interpolation_model(eb)
    initialize_melt_model(eb, case_inputs)
    initialize_topography_model(eb, case_inputs, ttype_grid)


def initialize_marker_output(eb: EarthBox) -> None:
    """ Initialize marker output.
    """
    model_manager = eb.model_manager
    marker_output = model_manager.config.output.marker_output
    keys = model_manager.config.output.marker_output_keys
    marker_output[keys.marker_age] = True
    marker_output[keys.marker_rho] = True
    marker_output[keys.marker_serpentinization] = True
    marker_output[keys.marker_extractable_meltfrac] = False
    marker_output[keys.marker_extracted_meltfrac] = False
    marker_output[keys.marker_strain_rate_plastic] = True


def initialize_staggered_grid(
        eb: EarthBox,
        ttype_grid: TtypeGrid
) -> None:
    """ Initialize staggered grid.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    staggered_grid = option_handler.staggered_grid
    option_names = staggered_grid.options.option_names
    key_names = staggered_grid.t_type_refinement_parameters.key_names
    option_handler.staggered_grid.initialize(
        pymodel,
        option_name=option_names.t_type_refined_grid,
        t_type_refinement_parameters={
            key_names.dx_highres: ttype_grid.dx_highres,
            key_names.dy_highres: ttype_grid.dy_highres,
            key_names.xo_highres: ttype_grid.xo_highres,
            key_names.xf_highres: ttype_grid.xf_highres,
            key_names.yf_highres: ttype_grid.yf_highres,
            }
            )


def initialize_boundary_conditions(
        eb: EarthBox,
        case_inputs: dict[str, float | int | bool]
) -> None:
    """ Initialize boundary conditions.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    celsius_to_kelvin = eb.units.celsius_to_kelvins
    cm_yr_to_m_s = eb.units.cm_per_yr_to_meters_per_second

    boundary_conditions = option_handler.boundary_conditions

    keys = CaseInputsKeys()

    option_names = boundary_conditions.options.option_names
    boundary_conditions.initialize(
        pymodel,
        option_name=option_names.lithospheric_extension_fixed_boundaries
        )
    boundary_conditions.pressure.initialize(pymodel, pressure_bc=1e5)

    (
        _temperature_base_lithosphere_warmer_initial_celsius,
        temperature_bottom_warmer_initial_celsius,
        temperature_bottom_cooler_transient_celsius
    ) = input_preprocessing.calculate_thermal_bcs_for_transient_cooling(case_inputs)

    boundary_conditions.temperature.initialize(
        pymodel,
        temperature_top=celsius_to_kelvin(0.0),
        temperature_bottom=celsius_to_kelvin(
            temperature_bottom_warmer_initial_celsius)
        )

    boundary_conditions.temperature.transient_bottom_temperature.initialize(
        pymodel,
        use_bottom_transient=case_inputs[keys.use_bottom_transient],
        temperature_bottom_transient_celsius=temperature_bottom_cooler_transient_celsius,
        start_time_bottom_transient_myr=case_inputs[keys.start_time_bottom_transient_myr],
        end_time_bottom_transient_myr=case_inputs[keys.end_time_bottom_transient_myr]
        )

    boundary_conditions.velocity.initialize(
        pymodel,
        full_velocity_extension_meters_per_second=cm_yr_to_m_s(
            case_inputs[keys.full_velocity_extension_cm_yr])
            )

    (
        time_step_parameters
    ) = input_preprocessing.get_time_step_parameters_for_velocity_stepping(case_inputs)
    boundary_conditions.velocity_step.initialize(
        pymodel,
        use_velocity_step=True,
        velocity_step_factor=time_step_parameters['velocity_step1_factor'],
        timestep_adjustment_factor=time_step_parameters['timestep_adjustment_factor_step1'],
        velocity_step_time=case_inputs[keys.velocity_step1_time_myr],
        velocity_second_step_factor=time_step_parameters['velocity_step2_factor'],
        timestep_second_adjustment_factor=time_step_parameters['timestep_adjustment_factor_step2'],
        velocity_second_step_time=case_inputs[keys.velocity_step2_time_myr]
        )


def initialize_marker_coordinates(eb: EarthBox) -> None:
    """ Initialize marker coordinates.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    markers = option_handler.markers

    marker_coordinates = markers.marker_coordinates
    option_names = marker_coordinates.options.option_names
    marker_coordinates.initialize(
        pymodel=pymodel, option_name=option_names.randomized)


def initialize_geometry(eb: EarthBox, case_inputs) -> None:
    """ Initialize geometry.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    keys = CaseInputsKeys()
    option_handler.geometry.sticky_air_geometry.initialize(
        pymodel, thick_air_meters=case_inputs[keys.thickness_sticky_meters])

    thicknesses = calculate_layer_thicknesses_meters(case_inputs)

    option_handler.geometry.earth_layering.initialize(
        pymodel,
        thick_upper_crust_meters=thicknesses['upper_continental_crust'],
        thick_lower_crust_meters=thicknesses['lower_continental_crust'],
        thick_upper_lith_meters=thicknesses['mantle_continental_lithosphere']/3.0,
        thick_middle_lith_meters=thicknesses['mantle_continental_lithosphere']/3.0,
        thick_lower_lith_meters=thicknesses['mantle_continental_lithosphere']/3.0,
        )

    option_handler.geometry.litho_strong_zones.initialize(
        pymodel,
        x_left_strong_meters=case_inputs[keys.xo_highres],
        x_right_strong_meters=case_inputs[keys.xf_highres],
        )


def calculate_layer_thicknesses_meters(
        case_inputs: dict[str, float | int]
) -> dict[str, float]:
    """ Calculate layer thicknesses.
    """
    keys = CaseInputsKeys()
    thickness_lithosphere_meters = case_inputs[
        keys.thickness_lithosphere_meters]
    thickness_continental_crustal_meters = case_inputs[
        keys.thickness_continental_crust_meters]
    thickness_upper_continental_crust_meters = case_inputs[
        keys.thickness_upper_continental_crust_meters]
    thickness_lower_continental_crust_meters = (
        thickness_continental_crustal_meters
        - thickness_upper_continental_crust_meters
        )
    thickness_mantle_continental_lithosphere_meters = (
        thickness_lithosphere_meters
        - thickness_continental_crustal_meters
        )
    return {
        'upper_continental_crust':
        thickness_upper_continental_crust_meters,
        'lower_continental_crust':
        thickness_lower_continental_crust_meters,
        'mantle_continental_lithosphere':
        thickness_mantle_continental_lithosphere_meters
        }


def initialize_marker_materials(
        eb: EarthBox,
        case_inputs: dict[str, float | int]
) -> None:
    """ Initialize marker materials.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    marker_materials = option_handler.markers.marker_materials

    materials_input_dict, lib_path = get_materials_input_dict()

    keys = CaseInputsKeys()

    option_names = marker_materials.options.option_names
    marker_materials.initialize(
        pymodel=pymodel,
        option_name=option_names.lithospheric_extension_lateral_strong_zones,
        paths={'materials_library_file': lib_path},
        materials_input_dict=materials_input_dict,
        viscosity_minimum_pascal_seconds=1e18,
        viscosity_maximum_pascal_seconds=1e26,
        use_fluid_pressure_for_yield=True,
        plastic_healing_rate=case_inputs[keys.plastic_healing_rate]
        )

    marker_materials.stress_limits.initialize(
        pymodel,
        powerlaw_stress_minimum_pascals=1e4,
        yield_stress_minimum_pascals=1e6,
        yield_stress_maximum_pascals=1e32
        )

    marker_materials.viscous_strain_softening_model.initialize(
        pymodel,
        use_viscous_strain_softening=False,
        viscous_strain_softening_factor=30.0
        )

    eb.model_manager.material_override.override_material_properties(**case_inputs)


def initialize_marker_temperature(
        eb: EarthBox,
        case_inputs: dict[str, float | int]
) -> None:
    """ Initialize marker temperature.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel
    celsius_to_kelvin = eb.units.celsius_to_kelvins

    markers = option_handler.markers
    marker_temperature = markers.marker_temperature

    option_names = marker_temperature.options.option_names
    key_names = marker_temperature.analytical_three_layer_parameters.key_names

    keys = CaseInputsKeys()

    (
        temperature_base_lithosphere_warmer_initial_celsius,
        _temperature_bottom_warmer_initial_celsius,
        _temperature_bottom_cooler_transient_celsius
    ) = input_preprocessing.calculate_thermal_bcs_for_transient_cooling(case_inputs)

    marker_temperature.initialize(
        pymodel=pymodel,
        option_name=option_names.analytical_three_layer,
        parameters_analytical_three_layer={
            key_names.temperature_base_lithosphere_kelvins: (
                celsius_to_kelvin(temperature_base_lithosphere_warmer_initial_celsius)
                ),
            key_names.amplitude_perturbation_meters: (
                case_inputs[keys.amplitude_perturbation_meters]),
            key_names.width_perturbation_meters: 10_000.0,
            key_names.conductivity_upper_crust_W_m_K: 2.25,
            key_names.conductivity_lower_crust_W_m_K: 2.0,
            key_names.conductivity_mantle_W_m_K: 2.0,
            key_names.heat_production_upper_crust_W_m3: (
                case_inputs[keys.radiogenic_heat_production_felsic_crust]),
            key_names.heat_production_lower_crust_W_m3: (
                case_inputs[keys.radiogenic_heat_production_mafic_crust]),
            key_names.heat_production_mantle_W_m3: 0.0,
            key_names.thick_thermal_lithosphere_meters: (
                case_inputs[keys.thickness_lithosphere_meters]),
            key_names.adiabatic_gradient: case_inputs[keys.adiabatic_temperature_gradient]
            }
        )


def initialize_marker_plasticity(eb: EarthBox) -> None:
    """ Initialize marker plasticity.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    markers = option_handler.markers

    option_names = markers.marker_friction.options.option_names
    option_handler.markers.marker_friction.initialize(
        pymodel,
        option_name=option_names.regular,
        update_randomization_at_each_time_step=True,
        time_randomization_factor=10.0
        )

    option_handler.markers.cohesion_initialize(pymodel)
    option_handler.markers.preexponential_initialize(pymodel)


def initialize_rock_properties(
        eb: EarthBox,
        case_inputs: dict[str, float | int]
) -> None:
    """ Initialize rock properties.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    rock_props = option_handler.rock_props

    density_names = rock_props.density.options.option_names
    conductivity_names = rock_props.thermal_conductivity.options.option_names
    rhocp_names = rock_props.rhocp.options.option_names
    rock_props.initialize(
        pymodel,
        option_name_density=density_names.liao14,
        option_name_thermal_conductivity=conductivity_names.sekiguchi_waples,
        option_name_rhocp=rhocp_names.temperature_dependent_waples,
        use_melt_lens=False,
        use_sediment_porosity=True
        )

    keys = CaseInputsKeys()

    rock_props.hydrothermal_circulation_model.initialize(
        pymodel,
        use_hydrothermal=True,
        hydrothermal_smoothing_factor=0.75,
        hydrothermal_max_depth_meters=case_inputs[
            keys.hydrothermal_max_depth_meters],
        hydrothermal_max_temperature_celsius=600.0,
        hydrothermal_nusselt_number=case_inputs[
            keys.hydrothermal_nusselt_number],
        use_plastic_strain_rate_for_hydrothermal=case_inputs[
            keys.use_plastic_strain_rate_for_hydrothermal],
        hydrothermal_plastic_strain_rate_reference=case_inputs[
            keys.hydrothermal_plastic_strain_rate_reference],
        use_plastic_strain_for_hydrothermal=case_inputs[
            keys.use_plastic_strain_for_hydrothermal],
        hydrothermal_plastic_strain_reference=case_inputs[
            keys.hydrothermal_plastic_strain_reference],
        hydrothermal_decay_length=case_inputs[
            keys.hydrothermal_decay_length],
        hydrothermal_buffer_distance=case_inputs[
            keys.hydrothermal_buffer_distance],
        sediment_thickness_threshold=case_inputs[
            keys.sediment_thickness_threshold]
        )

    rock_props.serpentinization_model.initialize(
        pymodel,
        use_serpentinization=True,
        serpentinization_temperature_celsius=340.5,
        maximum_serpentinization_depth_meters=case_inputs[
            keys.maximum_serpentinization_depth_meters],
        maximum_serpentinization_rate=case_inputs[
            keys.maximum_serpentinization_rate],
        nominal_strain_rate_serpentinization=1e-13
        )


def initialize_stokes_continuity_solver(
        eb: EarthBox,
        case_inputs: dict[str, float | int]
) -> None:
    """ Initialize Stokes continuity solver.
    """
    keys = CaseInputsKeys()
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    sc_solver = option_handler.stokes_continuity_solver
    velocity_type_names = sc_solver.velocity_type.options.option_names

    sc_solver.initialize(
        pymodel,
        option_name=velocity_type_names.velocity_from_stokes_solver,
        gravity_x_meters_per_second_per_second=0.0,
        gravity_y_meters_per_second_per_second=9.8,
        use_interface_stabilization=True
        )

    global_plasticity_loop = option_handler.global_plasticity_loop
    global_loop_names = global_plasticity_loop.options.option_names

    global_plasticity_loop.initialize(
        pymodel,
        option_name=global_loop_names.nodal_plasticity_loop,
        tolerance_picard=case_inputs[keys.tolerance_picard],
        maximum_number_of_global_plasticity_iterations=(
            case_inputs[keys.maximum_number_of_global_plasticity_iterations])
        )


def initialize_heat_solver(eb: EarthBox) -> None:
    """ Initialize heat solver.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    option_handler.heat_solver.initialize(
        pymodel,
        solve_heat_equation=True,
        maximum_temperature_change_kelvins=70.0,
        use_adiabatic_heating=True,
        use_shear_heating=True,
        use_sticky_correction=True
        )


def initialize_advection_model(
        eb: EarthBox,
        case_inputs: dict[str, float | int]
) -> None:
    """ Initialize advection model.
    """
    keys = CaseInputsKeys()
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    advection = option_handler.advection
    option_names = advection.options.option_names
    advection.initialize(
        pymodel,
        option_name=option_names.runge_kutta_4th_order,
        marker_cell_displ_max=case_inputs[
            keys.marker_cell_displ_max],
        subgrid_diff_coef_temp=1.0,
        subgrid_diff_coef_stress=1.0
        )


def initialize_interpolation_model(eb: EarthBox) -> None:
    """ Initialize interpolation model.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    option_handler.interpolate.initialize(
        pymodel,
        interpolate_normal_viscosity_with_harmonic_avg=True,
        interpolate_temps_back_to_markers_at_startup=True
        )


def initialize_melt_model(
        eb: EarthBox,
        case_inputs: dict[str, float | int]
) -> None:
    """ Initialize melt model.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    melt_model = option_handler.melt_model

    keys = CaseInputsKeys()

    melt_model.initialize(
        pymodel,
        use_melting=True,
        use_melt_viscosity=True,
        use_melt_thermal_properties=True,
        viscosity_melt_pascal_seconds=1e18,
        use_exponential_viscosity_reduction=True,
        use_depletion_density=True,
        debug=False
        )

    melt_model.melt_damage_model.initialize(
        pymodel,
        use_melt_damage=case_inputs[keys.use_melt_damage],
        melt_damage_distance=case_inputs[keys.melt_damage_distance],
        melt_damage_factor=case_inputs[keys.melt_damage_factor],
        use_probabilistic_melt_damage=case_inputs[keys.use_probabilistic_melt_damage],
        maximum_damage_probability=case_inputs[keys.maximum_damage_probability],
        intermediate_damage_probability=case_inputs[keys.intermediate_damage_probability],
        magmatic_crust_height_threshold=case_inputs[keys.magmatic_crust_height_threshold],
        magmatic_crust_height_minimum=case_inputs[keys.magmatic_crust_height_minimum],
        magmatic_crust_height_intermediate=case_inputs[keys.magmatic_crust_height_intermediate],
        magmatic_crust_height_maximum=case_inputs[keys.magmatic_crust_height_maximum],
        )

    melt_model.extraction_model.initialize(
        pymodel,
        use_extraction=True,
        use_shallow_mantle_injection=True,
        use_random_injection_subdomain=True,
        use_normal_injection_subdomain=True,
        smoothing_radius_drainage=10_000.0,
        smoothing_radius_fractionation=10_000.0,
        characteristic_injection_width=10_000.0,
        mantle_search_width=200_000.0,
        number_of_injection_subdomains=10,
        magma_height_limit=1000.0,
        emplacement_temperature_celsius=1200.0,
        use_gabbroic_fractionation=True,
        fractionation_threshold_limit=2000.0,
        maximum_shallow_injection_depth=25_000.0,
        extraction_fraction=0.99
        )

    melt_model.extrusion_model.initialize(
        pymodel,
        use_extrusion=True,
        extrusion_volume_factor=case_inputs[
            keys.extrusion_volume_factor],
        extrusion_volume_factor_max=case_inputs[
            keys.extrusion_volume_factor_max],
        characteristic_magmatic_crust_height_min=case_inputs[
            keys.characteristic_magmatic_crust_height_min],
        characteristic_magmatic_crust_height_max=case_inputs[
            keys.characteristic_magmatic_crust_height_max],
        width_eruption_domain_fixed=case_inputs[
            keys.width_eruption_domain_fixed],
        width_eruption_domain_fixed_max=case_inputs[
            keys.width_eruption_domain_fixed_max],
        characteristic_flow_length_subaerial=case_inputs[
            keys.characteristic_flow_length_subaerial],
        characteristic_flow_length_submarine=case_inputs[
            keys.characteristic_flow_length_submarine],
        residual_lava_thickness_subaerial=case_inputs[
            keys.residual_lava_thickness_subaerial],
        residual_lava_thickness_submarine=case_inputs[
            keys.residual_lava_thickness_submarine],
        decimation_factor=4,
        use_random_eruption_location=case_inputs[
            keys.use_normal_eruption_location],
        use_normal_eruption_location=case_inputs[
            keys.use_normal_eruption_location],
        porosity_initial_lava_flow=case_inputs[
            keys.porosity_initial_lava_flow],
        decay_depth_lava_flow=2000.0,
        use_eruption_interval=case_inputs[
            keys.use_eruption_interval],
        eruption_interval_yr=case_inputs[
            keys.eruption_interval_yr],
        )


def initialize_topography_model(
        eb: EarthBox,
        case_inputs: dict[str | float | int],
        ttype_grid: TtypeGrid
) -> None:
    """ Initialize topography model.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    topography = option_handler.topography

    keys = CaseInputsKeys()

    advection_types = topography.topo_type.options.option_names
    topography.initialize(
        pymodel,
        use_topography=True,
        number_of_topography_nodes=ttype_grid.calculate_number_of_topo_nodes(),
        width_of_topo_grid_meters=ttype_grid.model_width_meters,
        option_name=advection_types.runge_kutta_with_interp,
        nsmooth_top_bottom=2,
        marker_search_factor=2.0
        )

    sealevel_types = topography.sealevel.options.option_names
    topography.sealevel.initialize(
        pymodel,
        option_name=sealevel_types.average_pressure,
        depth_of_sealevel_meters=case_inputs[keys.thickness_sticky_meters],
        base_level_shift_end_time_myr=case_inputs[keys.base_level_shift_end_time_myr],
        base_level_shift_meters=case_inputs[keys.base_level_shift_meters]
        )

    topography.sealevel.reference_lithosphere.initialize(
        pymodel,
        thickness_upper_continental_crust_meters=22_000.0,
        thickness_lower_continental_crust_meters=10_000.0,
        thickness_lithosphere_meters=125_000.0,
        gridy_spacing_meters=100.0,
        temperature_top_celsius=0.0,
        temperature_moho_celsius=600.0,
        temperature_base_lith_celsius=case_inputs[keys.temperature_base_lithosphere_celsius],
        adiabatic_gradient_kelvin_km=case_inputs[keys.adiabatic_temperature_gradient],
        )

    m2_yr_to_m2_s = eb.units.m2_yr_to_m2_s
    m_yr_to_m_s = eb.units.m_yr_to_m_s
    yr_to_s = eb.units.years_to_seconds
    mm_yr_to_m_s = eb.units.mm_per_yr_to_meters_per_seconds

    (
        sediment_transport_timestep_yr
    ) = input_preprocessing.calculate_sediment_transport_timestep(
        case_inputs[keys.full_velocity_extension_cm_yr],
        case_inputs[keys.dx_highres],
        case_inputs[keys.dy_highres],
        case_inputs[keys.marker_cell_displ_max],
        number_of_transport_timesteps_per_model_timestep=5
        )
    print(">> sediment_transport_timestep_yr:", sediment_transport_timestep_yr)

    topography.sediment_transport_model.initialize(
        pymodel,
        use_downhill_diffusion=True,
        subaerial_slope_diffusivity=m2_yr_to_m2_s(
            case_inputs[keys.subaerial_slope_diffusivity]),
        precipitation_rate=m_yr_to_m_s(case_inputs[keys.precipitation_rate]),
        subaerial_transport_coefficient=(
            case_inputs[keys.subaerial_transport_coefficient]),
        submarine_slope_diffusivity=m2_yr_to_m2_s(
            case_inputs[keys.submarine_slope_diffusivity]),
        submarine_diffusion_decay_depth=case_inputs[
            keys.submarine_diffusion_decay_depth],
        transport_timestep=yr_to_s(sediment_transport_timestep_yr),
        pelagic_sedimentation_rate_meters_per_second=mm_yr_to_m_s(
            case_inputs[keys.pelagic_sedimentation_rate]),
        use_compaction_correction=True
        )


def activate_loop_plots(eb: EarthBox) -> None:
    """ Activate loop plots.
    """
    eb.model_manager.model_plots_2d.activate_scalar_loop_plot(
        scalar_name=eb.model_manager.model_plots_2d.scalar_names.temperature,
        contour_interval=100.0,
        minimum_value=0.0,
        maximum_value=1400.0,
        excluded_value=-1e38,
        plot_contours=True,
        grid_plot_type='nomesh'
        )

    plot_types = MarkerPlotNames
    eb.model_manager.model_plots_2d.activate_marker_loop_plot(
        plot_type=plot_types.Composition.name,
        plot_topography=True,
        plot_mesh=False,
        marker_size=0.25
        )


if __name__ == '__main__':
    run_model(use_loop_plotting=False, model_case_name=get_model_case_name())
