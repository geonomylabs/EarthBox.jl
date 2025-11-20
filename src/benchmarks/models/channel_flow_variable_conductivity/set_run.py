#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Script for working with models where parameters are input through the API.

Usage
-----
To run model:
    % python set_run.py run_model
To plot markers:
    % python set_run.py plot_markers
To plot scalars:
    % python set_run.py plot_scalars
To plot velocity:
    % python set_run.py plot_velocity

Note that plots can be made in time loops by activating loop plotting using
the setup_loop_plotting function.

"""
from typing import Callable
import os
import sys
from earthbox import EarthBox
from earthbox.plot_tools import ModelPlots2D
from earthbox.plot_tools import ScalarNames
from earthbox.plot_tools import MarkerPlotNames
from earthbox import MaterialsLibrary


def get_materials_lib_path():
    """ Get materials library path.
    """
    lib = MaterialsLibrary()
    return lib.benchmarks.path


def get_output_path():
    """ Get output path.
    """
    return r'../../results/channel_flow_variable_conductivity_output'


def get_materials_input_dict() -> dict[str, dict[str, float | int | str]]:
    """ Get materials input dict.
    """
    lib = MaterialsLibrary()
    mat_names = lib.benchmarks.materials
    materials_input_dict = {
        0: {
            'mat_name':mat_names.
            constant_viscosity_channel_variable_thermal_conductivity,
            'mat_type': 'General',
            'mat_domain': 'GeneralDomain',
            'red_fraction': 0.0,
            'green_fraction': 0.0,
            'blue_fraction': 1.0,
            }
        }
    return materials_input_dict


def setup_loop_plotting(eb: EarthBox, output_dir: str) -> None:
    """ Setup plotting for time loops.
    """
    # Initialize plotting parameters for model initialization and time loop
    eb.model_manager.model_plots_2d.initialize(
            plot_output_path=os.path.join(output_dir, 'plots'),
            material_library_file_path=get_materials_lib_path(),
            materials_input_dict=get_materials_input_dict(),
            iplot_contour_labels=1,
            color_map='bwr',
            figure_dpi=150.0,
            plot_length_units='km',
            plot_time_units='Myr',
            plot_velocity_units='cm/yr',
            axis_fraction_for_color_bar=0.06,
            axis_fraction_for_color_bar_gap=0.12,
            dimensions=(0.0, 20.0, 0.0, 12.5),
            use_close_up=False,
            dim_close_up=(-5.0, 95.0, 250.0, 305.0),
            xyspacing=(2.0, 1.0)
            )
    # Activate loop plot for temperature and define plot parameters
    eb.model_manager.model_plots_2d.activate_scalar_loop_plot(
        scalar_name=eb.model_manager.model_plots_2d.scalar_names.temperature,
        contour_interval=100.0,
        minimum_value=0.0,
        maximum_value=1300.0,
        excluded_value=-1e38,
        plot_contours=True,
        grid_plot_type='nomesh'
        )
    # Activate loop plot for markers and define plot parameters
    eb.model_manager.model_plots_2d.activate_marker_loop_plot(
        plot_type='Composition',
        plot_topography=False,
        plot_mesh=False,
        marker_size=0.25
        )


def setup_model() -> EarthBox:
    """ Setup model
    """
    # Define the output directory
    output_dir = get_output_path()
    # Initialize EarthBox model
    eb = EarthBox(
        paths={'output_dir': output_dir},
        number_of_basic_nodes_x=51,
        number_of_basic_nodes_y=11,
        model_width_meters=3e4,
        model_height_meters=1.25e4,
        number_of_markers_per_cell_x=5.0,
        number_of_markers_per_cell_y=5.0
        )
    # Get conversion functions
    celsius_to_kelvin = eb.units.celsius_to_kelvins
    # cm_yr_to_m_s = eb.units.cm_per_yr_to_meters_per_second
    # Get the model manager object. The model manager provides access to the
    # option handler and pymodel.
    model_manager = eb.model_manager
    # The option handler provides access to different components of the model
    # and there associated options.
    option_handler = model_manager.option_handler
    # The PyModelData class is a wrapper for a Numba structref data container
    # that provides access to all model parameters and arrays. Model data is
    # organized into different groups and for each group there is a
    # corresponding data class for parameters and arrays. Parameter objects
    # have attributes for name, units, value and description. Array objects
    # have attributes for name, units, array and description.
    pymodel = model_manager.pymodel
    # Initialize the basic grid using permitted option names
    option_names = option_handler.staggered_grid.options.option_names
    option_handler.staggered_grid.initialize(
        pymodel, grid_type=option_names.uniform_grid
        )
    # Initialize the boundary conditions using permitted option names
    option_names = option_handler.boundary_conditions.options.option_names
    option_handler.boundary_conditions.initialize(
        pymodel,
        model_type=option_names.vertical_channel_flow_shear_heating,
        pressure_pascals=3e7,
        temperature_top_kelvins=celsius_to_kelvin(25.0),
        temperature_bottom_kelvins=celsius_to_kelvin(25.0)
        )

    markers = option_handler.markers
    # Initialize the marker coordinates using permitted option names
    marker_coordinates = markers.marker_coordinates
    option_names = marker_coordinates.options.option_names
    marker_coordinates.initialize(
        pymodel=pymodel, marker_distribution=option_names.regular)
    # Initialize the marker materials using permitted option names
    marker_materials = markers.marker_materials
    option_names = marker_materials.options.option_names
    marker_materials.initialize(
        pymodel=pymodel,
        material_model=option_names.uniform_composition,
        paths={'materials_library_file': get_materials_lib_path()},
        materials_input_dict=get_materials_input_dict(),
        viscosity_minimum_pascal_seconds=1e18,
        viscosity_maximum_pascal_seconds=1e25,
        )
    # Initialize the marker temperature using permitted option names.
    marker_temperature = markers.marker_temperature
    option_names = marker_temperature.options.option_names
    key_names = marker_temperature.linear_parameters.key_names
    marker_temperature.initialize(
        pymodel=pymodel,
        initial_temperature=option_names.linear,
        parameters_linear={
            key_names.temperature_top_kelvins: celsius_to_kelvin(25.0),
            key_names.temperature_bottom_kelvins: celsius_to_kelvin(25.0)
            }
        )

    # Initialize marker density and thermal conductivity models using permitted
    # option names
    rock_props = option_handler.rock_props
    density_names = rock_props.density.options.option_names
    cond_names = rock_props.thermal_conductivity.options.option_names
    rock_props.initialize(
        pymodel,
        density_model=density_names.liao14,
        thermal_conductivity_model=cond_names.temperature_dependent,
        )

    # Initialize the Stokes-continuity solver (sc_solver) using permitted
    # option names
    sc_solver = option_handler.stokes_continuity_solver
    global_loop_names = sc_solver.global_plasticity_loop.options.option_names
    velocity_type_names = sc_solver.velocity_type.options.option_names
    marker_plasticity_names = sc_solver.marker_plasticity.options.option_names
    sc_solver.initialize(
        pymodel,
        velocity_type=velocity_type_names.velocity_from_stokes_solver,
        global_plasticity_loop=global_loop_names.nodal_plasticity_loop,
        tolerance_picard=0.001,
        maximum_number_of_global_plasticity_iterations=1,
        marker_plasticity_model=
        marker_plasticity_names.pure_elastic_stress_forecast,
        gravity_x_meters_per_second_per_second=0.0,
        gravity_y_meters_per_second_per_second=0.0
        )

    # Initialize the heat solver
    option_handler.heat_solver.initialize(
        pymodel,
        solve_heat_equation=True,
        maximum_temperature_change_kelvins=10.0,
        use_adiabatic_heating=False,
        use_shear_heating=True,
        )

    # Initialize the marker advection method using permitted option names
    advection = option_handler.advection
    option_names = advection.options.option_names
    advection.initialize(
        pymodel,
        advection_scheme=option_names.no_motion,
        maximum_marker_cell_displacement_fraction=50.0,
        subgrid_temperature_diffusion_coefficient=1.0,
        subgrid_stress_diffusion_coefficient=1.0
        )

    # Initialize the marker interpolation method
    option_handler.interpolate.initialize(
        pymodel,
        interpolate_normal_viscosity_with_harmonic_avg=False,
        interpolate_temps_back_to_markers_at_startup=True
        )

    return eb


def set_and_run_model() -> None:
    """ Set up model and run model time steps
    """
    # Set up model using the Earthbox API
    eb = setup_model()
    # Configure solver options
    eb.model_manager.config.solver.use_mumps = False
    eb.model_manager.config.solver.nproc = 1
    # Get conversion function
    years_to_seconds = eb.units.years_to_seconds
    # Run model time steps
    eb.run_time_steps(
        number_of_time_steps=500,
        viscoelastic_time_step_seconds=years_to_seconds(1e5),
        output_time_step_seconds=years_to_seconds(5e5),
        make_backup=False
        )


def get_model_plots_2d():
    """ Get model plots 2D.
    """
    lib = MaterialsLibrary()
    mp2d = ModelPlots2D(
        plot_output_path=os.path.join(get_output_path(), 'plots'),
        material_library_file_path=lib.benchmarks.path,
        materials_input_dict=get_materials_input_dict(),
        model_output_path=get_output_path(),
        nsteps=30,
        nskip=0,
        istart=0,
        iplot_contour_labels=1,
        color_map='bwr',
        figure_dpi=150.0,
        plot_length_units='km',
        plot_time_units='Myr',
        plot_velocity_units='cm/yr',
        axis_fraction_for_color_bar=0.06,
        axis_fraction_for_color_bar_gap=0.12,
        dimensions=(0.0, 30.0, 0.0, 12.5),
        use_close_up=False,
        dim_close_up=(-5.0, 95.0, 250.0, 305.0),
        xyspacing=(2.0, 2.0)
        )
    return mp2d


def plot_scalars():
    """ Plot scalars.
    """
    mp2d = get_model_plots_2d()
    scalar_names = ScalarNames()
    mp2d.plot_scalars(
        scalar_names.temperature,
        contour_interval=100.0,
        minimum_value=0.0,
        maximum_value=1300.0,
        plot_contours=True,
        grid_plot_type='nomesh'
        )
    mp2d.plot_scalars(
        scalar_names.viscosity,
        contour_interval=0.5,
        minimum_value=19.0,
        maximum_value=23.0,
        plot_contours=True,
        grid_plot_type='nomesh'
        )

def plot_velocity():
    """ Plot velocity.
    """
    mp2d = get_model_plots_2d()
    mp2d.plot_velocity()


def plot_markers():
    """ Plot markers.
    """
    mp2d = get_model_plots_2d()
    plot_types = MarkerPlotNames()
    mp2d.plot_markers(
        plot_type=plot_types.composition,
        plot_topography=False,
        plot_mesh=False,
        marker_size=1.0
        )


def get_func(option_name: str) -> Callable[[None], None]:
    """ Get function based on option name.
    """
    func_dict = {
        'run_model': set_and_run_model,
        'plot_markers': plot_markers,
        'plot_scalars': plot_scalars,
        'plot_velocity': plot_velocity
        }
    if option_name not in func_dict:
        raise ValueError(
            f'{option_name} is not a valid option. '
            f'Valid option are {list(func_dict.keys())}'
            )
    return func_dict[option_name]


if __name__ == '__main__':
    get_func(sys.argv[1])()
