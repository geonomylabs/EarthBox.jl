#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Script for working with models where parameters are input through the API.

Usage
-----
To run model:
    % python run.py

Note that plots can be made in time loops by activating loop plotting using
the setup_loop_plotting function.

"""
from paths import get_output_path
from materials import get_materials_input_dict
from plot import get_plot_parameters
from earthbox import EarthBox
from earthbox.plot_tools.utils.loop_plots import setup_loop_plotting


def run_model(use_loop_plotting: bool=False) -> None:
    """ Set up model and run model time steps
    """
    output_dir = get_output_path()
    eb = setup_model(output_dir)

    eb.model_manager.config.solver.use_mumps = True
    eb.model_manager.config.solver.nproc = 10
    eb.model_manager.config.solver.pymumps_timeout = 60.0

    if use_loop_plotting:
        setup_loop_plotting(eb, activate_loop_plots, get_plot_parameters)

    years_to_seconds = eb.units.years_to_seconds
    eb.run_time_steps(
        number_of_time_steps=2500,
        viscoelastic_time_step_seconds=years_to_seconds(50_000),
        output_time_step_seconds=years_to_seconds(500_000),
        make_backup=True
        )


def setup_model(output_dir: str) -> EarthBox:
    """ Setup model
    """

    # average x-resolution = 1_00_000/1850 = 540 m
    # average y-resolution = 160_000/400 = 320 m
    eb = EarthBox(
        restart_from_backup=False,
        paths={'output_dir': output_dir},
        xnum=1851,
        ynum=401,
        xsize=1_000_000.0,
        ysize=160_000.0,
        nmarkers_cell_x=8.0,
        nmarkers_cell_y=5.0
        )
    initialize_model_input(eb)
    return eb


def initialize_model_input(eb: EarthBox) -> None:
    """ User defined functions to initialize input via the EarthBox API.
    """
    initialize_marker_output(eb)
    initialize_staggered_grid(eb)
    initialize_boundary_conditions(eb)
    initialize_marker_coordinates(eb)
    initialize_geometry(eb)
    initialize_marker_materials(eb)
    initialize_marker_temperature(eb)
    initialize_marker_plasticity(eb)
    initialize_rock_properties(eb)
    initialize_stokes_continuity_solver(eb)
    initialize_heat_solver(eb)
    initialize_advection_model(eb)
    initialize_interpolation_model(eb)
    initialize_melt_model(eb)
    initialize_topography_model(eb)


def initialize_marker_output(eb: EarthBox) -> None:
    """ Initialize marker output.
    """
    model_manager = eb.model_manager
    marker_output = model_manager.config.output.marker_output
    marker_output_keys = model_manager.config.output.marker_output_keys
    marker_output[marker_output_keys.marker_extractable_meltfrac] = False
    marker_output[marker_output_keys.marker_extracted_meltfrac] = False
    marker_output[marker_output_keys.marker_serpentinization] = True
    marker_output[marker_output_keys.marker_strain_rate_plastic] = False
    marker_output[marker_output_keys.marker_TK] = True
    marker_output[marker_output_keys.marker_age] = True


def initialize_staggered_grid(eb: EarthBox) -> None:
    """ Initialize staggered grid.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    staggered_grid = option_handler.staggered_grid
    option_names = staggered_grid.options.option_names
    key_names = staggered_grid.t_type_refinement_parameters.key_names

    # number of basic cells in high-res domain in x-direction = 300_000/200 = 1500
    # number of basic cells in high-res domain in y-direction = 70_000/200 = 350
    option_handler.staggered_grid.initialize(
        pymodel, option_name=option_names.t_type_refined_grid,
        t_type_refinement_parameters={
            key_names.dx_highres: 200.0,
            key_names.dy_highres: 200.0,
            key_names.xo_highres: 350_000.0,
            key_names.xf_highres: 650_000.0,
            key_names.yf_highres: 70_000.0,
            }
        )


def initialize_boundary_conditions(eb: EarthBox) -> None:
    """ Initialize boundary conditions.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    celsius_to_kelvin = eb.units.celsius_to_kelvins
    cm_yr_to_m_s = eb.units.cm_per_yr_to_meters_per_second

    boundary_conditions = option_handler.boundary_conditions

    option_names = boundary_conditions.options.option_names
    boundary_conditions.initialize(
        pymodel,
        option_name=option_names.lithospheric_extension_fixed_boundaries
        )

    boundary_conditions.pressure.initialize(pymodel, pressure_bc=1e5)

    boundary_conditions.temperature.initialize(
        pymodel,
        temperature_top=celsius_to_kelvin(0.0),
        temperature_bottom=celsius_to_kelvin(1330.0 + 0.4*25.0)
        )

    boundary_conditions.velocity.initialize(
        pymodel,
        full_velocity_extension_meters_per_second=cm_yr_to_m_s(0.2)
        )

    # Minimum resolution is 200 m. At 2 cm/yr maximum time step should be
    # 200/0.02 = 10_000 years.
    boundary_conditions.velocity_step.initialize(
        pymodel,
        use_velocity_step=True,
        velocity_step_factor=1.5/0.2, # 0.2 cm/yr to 1.5 cm/yr
        timestep_adjustment_factor=13_000.0/50_000.0, # 50_000.0 to 13_000.0 yr
        velocity_step_time=40.0, # Myr
        velocity_second_step_factor=2.0/1.5, # 1.5 cm/yr to 2 cm/yr
        timestep_second_adjustment_factor=10_000.0/13_000.0, # 13_000 yr to 10_000 yr
        velocity_second_step_time=50.0
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


def initialize_geometry(eb: EarthBox) -> None:
    """ Initialize geometry.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    option_handler.geometry.sticky_air_geometry.initialize(
        pymodel, thick_air_meters=10_000.0)

    thicknesses = calculate_layer_thicknesses_meters()

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
        x_left_strong_meters=350_000.0,
        x_right_strong_meters=650_000.0
        )

def calculate_layer_thicknesses_meters():
    """ Calculate layer thicknesses.
    """
    thickness_lithosphere_meters = 125_000.0
    thickness_continental_crustal_meters = 35_000.0
    thickness_upper_continental_crust_meters = 22_000.0
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

def initialize_marker_materials(eb: EarthBox) -> None:
    """ Initialize marker materials.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    marker_materials = option_handler.markers.marker_materials

    materials_input_dict, lib_path = get_materials_input_dict()

    option_names = marker_materials.options.option_names
    marker_materials.initialize(
        pymodel=pymodel,
        option_name=option_names.lithospheric_extension_lateral_strong_zones,
        paths={'materials_library_file': lib_path},
        materials_input_dict=materials_input_dict,
        viscosity_minimum_pascal_seconds=1e18,
        viscosity_maximum_pascal_seconds=1e26
        )

    marker_materials.stress_limits.initialize(
        pymodel,
        powerlaw_stress_minimum_pascals=1e4,
        yield_stress_minimum_pascals=0.0,
        yield_stress_maximum_pascals=1e32
        )

    marker_materials.viscous_strain_softening_model.initialize(
        pymodel,
        use_viscous_strain_softening=False,
        viscous_strain_softening_factor=30.0
        )


def initialize_marker_temperature(eb: EarthBox) -> None:
    """ Initialize marker temperature.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel
    celsius_to_kelvin = eb.units.celsius_to_kelvins

    markers = option_handler.markers
    marker_temperature = markers.marker_temperature

    option_names = marker_temperature.options.option_names
    key_names = marker_temperature.analytical_three_layer_parameters.key_names

    marker_temperature.initialize(
        pymodel=pymodel,
        option_name=option_names.analytical_three_layer,
        parameters_analytical_three_layer={
            key_names.temperature_base_lithosphere_kelvins: celsius_to_kelvin(1330.0),
            key_names.amplitude_perturbation_meters: 5_000.0,
            key_names.width_perturbation_meters: 10_000.0,
            key_names.conductivity_upper_crust_W_m_K: 2.25,
            key_names.conductivity_lower_crust_W_m_K: 2.0,
            key_names.conductivity_mantle_W_m_K: 2.0,
            key_names.heat_production_upper_crust_W_m3: 1.8e-6,
            key_names.heat_production_lower_crust_W_m3: 0.5e-6,
            key_names.heat_production_mantle_W_m3: 0.0,
            key_names.thick_thermal_lithosphere_meters: 125_000.0,
            key_names.adiabatic_gradient: 0.4
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


def initialize_rock_properties(eb: EarthBox) -> None:
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
        use_sediment_porosity=True
        )

    rock_props.hydrothermal_circulation_model.initialize(
        pymodel,
        use_hydrothermal=True,
        hydrothermal_smoothing_factor=0.75,
        hydrothermal_nusselt_number=4.0,
        hydrothermal_max_temperature_celsius=600.0,
        hydrothermal_max_depth_meters=6000.0
        )

    rock_props.serpentinization_model.initialize(
        pymodel,
        use_serpentinization=True,
        serpentinization_temperature_celsius=340.5,
        maximum_serpentinization_depth_meters=20_000.0,
        maximum_serpentinization_rate=1e-11,
        nominal_strain_rate_serpentinization=1e-13
        )


def initialize_stokes_continuity_solver(eb: EarthBox) -> None:
    """ Initialize Stokes continuity solver.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    sc_solver = option_handler.stokes_continuity_solver
    global_loop_names = sc_solver.global_plasticity_loop.options.option_names
    velocity_type_names = sc_solver.velocity_type.options.option_names
    marker_plasticity_names = sc_solver.marker_plasticity.options.option_names
    sc_solver.initialize(
        pymodel,
        option_name=velocity_type_names.velocity_from_stokes_solver,
        global_plasticity_loop=global_loop_names.nodal_plasticity_loop,
        tolerance_picard=0.01,
        maximum_number_of_global_plasticity_iterations=3,
        marker_plasticity_model=marker_plasticity_names.viscoelastic_stress_forecast,
        gravity_x_meters_per_second_per_second=0.0,
        gravity_y_meters_per_second_per_second=9.8,
        use_interface_stabilization=True
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


def initialize_advection_model(eb: EarthBox) -> None:
    """ Initialize advection model.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    advection = option_handler.advection
    option_names = advection.options.option_names
    advection.initialize(
        pymodel,
        option_name=option_names.runge_kutta_4th_order,
        marker_cell_displ_max=1.0,
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


def initialize_melt_model(eb: EarthBox) -> None:
    """ Initialize melt model.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    melt_model = option_handler.melt_model

    melt_model.initialize(
        pymodel,
        use_melting=True,
        use_melt_viscosity=True,
        use_melt_thermal_properties=True,
        viscosity_melt_pascal_seconds=1e18,
        use_depletion_density=True
        )

    melt_model.extraction_model.initialize(
        pymodel,
        use_extraction=True,
        use_shallow_mantle_injection=True,
        use_random_injection_subdomain=True,
        smoothing_radius_drainage=10_000.0,
        characteristic_injection_width=10_000.0,
        number_of_injection_subdomains=10,
        magma_height_limit=2500.0,
        emplacement_temperature_celsius=1200.0,
        use_gabbroic_fractionation=True,
        fractionation_threshold_limit=2000.0,
        maximum_shallow_injection_depth=40_000.0,
        extraction_fraction=0.99
        )

    melt_model.extrusion_model.initialize(
        pymodel,
        use_extrusion=True,
        extrusion_volume_factor=0.1,
        characteristic_flow_length_subaerial=100_000.0,
        characteristic_flow_length_submarine=2000.0,
        residual_lava_thickness_subaerial=30.0,
        residual_lava_thickness_submarine=10.0,
        decimation_factor=4,
        use_random_eruption_location=True,
        use_normal_eruption_location=True
        )


def initialize_topography_model(eb: EarthBox) -> None:
    """ Initialize topography model.
    """
    option_handler = eb.model_manager.option_handler
    pymodel = eb.model_manager.pymodel

    topography = option_handler.topography
    advection_types = topography.topo_type.options.option_names
    sealevel_types = topography.sealevel.options.option_names
    topography.initialize(
        pymodel,
        use_topography=True,
        number_of_topography_nodes=5001,
        width_of_topo_grid_meters=1_000_000.0,
        option_name=advection_types.runge_kutta_with_interp,
        nsmooth_top_bottom=2,
        marker_search_factor=2.0
        )

    topography.sealevel.initialize(
        pymodel,
        depth_of_sealevel_meters=10_000.0,
        option_name=sealevel_types.average_pressure,
        base_level_shift_meters=0.0,
        base_level_shift_end_time_myr=0.0
        )

    topography.sealevel.reference_lithosphere.initialize(
        pymodel,
        thickness_upper_continental_crust_meters=22_000.0,
        thickness_lower_continental_crust_meters=10_000.0,
        thickness_lithosphere_meters=125_000.0,
        gridy_spacing_meters=10.0,
        temperature_top_celsius=0.0,
        temperature_moho_celsius=600.0,
        temperature_base_lith_celsius=1330.0,
        adiabatic_gradient_kelvin_km=0.4,
        use_linear_segments=False
        )

    m2_yr_to_m2_s = eb.units.m2_yr_to_m2_s
    m_yr_to_m_s = eb.units.m_yr_to_m_s
    yr_to_s = eb.units.years_to_seconds
    mm_yr_to_m_s = eb.units.mm_per_yr_to_meters_per_seconds

    topography.sediment_transport_model.initialize(
        pymodel,
        use_downhill_diffusion=True,
        subaerial_slope_diffusivity=m2_yr_to_m2_s(0.25),
        precipitation_rate=m_yr_to_m_s(1.0),
        subaerial_transport_coefficient=1.0e-4,
        submarine_slope_diffusivity=m2_yr_to_m2_s(1.0),
        submarine_diffusion_decay_depth=1000.0,
        transport_timestep=yr_to_s(5_000.0),
        pelagic_sedimentation_rate_meters_per_second=mm_yr_to_m_s(0.0),
        use_compaction_correction=True
        )

    topography.salt_deposition_model.initialize(
        pymodel,
        use_salt_deposition=False,
        salt_start_time_myr=0.0,
        salt_end_time_myr=0.0,
        salt_deposition_rate_meters_per_second=0.0
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

    eb.model_manager.model_plots_2d.activate_marker_loop_plot(
        plot_type='Composition',
        plot_topography=True,
        plot_mesh=False,
        marker_size=0.25
        )


if __name__ == '__main__':
    run_model(use_loop_plotting=False)
