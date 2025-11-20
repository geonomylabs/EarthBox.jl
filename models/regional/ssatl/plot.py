""" Plotting functions.

Usage
-----
To plot markers:
    % python plot.py plot_markers model_case_name istart

To plot scalars:
    % python plot.py plot_scalars model_case_name istart

To plot velocity:
    % python plot.py plot_velocity model_case_name istart

"""
import os
import sys
from enum import Enum
from typing import Callable
from materials import get_materials_input_dict
from paths import get_output_path
from case_inputs import get_case_inputs
from earthbox.input_tools.case_inputs_tools import CaseInputsKeys
from earthbox.input_tools import input_preprocessing
from earthbox.input_tools.get_inputs import get_istart_iend
from earthbox.plot_tools import ModelPlots2D
from earthbox.plot_tools import ScalarNames
from earthbox.plot_tools import MarkerPlotNames
from earthbox.plot_tools.utils.color_maps import MatplotlibColorMapNames
from earthbox.plot_tools.log_file_plots.magmatic_crust import plot_magmatic_crust_height
from earthbox.utils.get_paths import get_model_case_name_for_plotting


class LocalPlotTypes(Enum):
    """ Local user defined plot types
    """
    GENERAL_ZOOM_OUT = 0
    ZOOM_IN_TO_HIGHRES = 2
    ZOOM_IN_2_HIGHRES_HF_GRAV = 3


def define_local_plot_type():
    """ Define local plot type

    The user should define the plot type here. This is used to define the
    specific plot parameters.

    """
    local_plot_type = LocalPlotTypes.ZOOM_IN_TO_HIGHRES.name
    #local_plot_type = LocalPlotTypes.ZOOM_IN_2_HIGHRES_HF_GRAV.name
    return local_plot_type


def get_marker_size(local_plot_type: str) -> float:
    """ Define zoom options and return selected option """
    marker_size_dict = {
        LocalPlotTypes.ZOOM_IN_TO_HIGHRES.name: 0.5,
        LocalPlotTypes.ZOOM_IN_2_HIGHRES_HF_GRAV.name: 0.25
        }
    return marker_size_dict[local_plot_type]


def get_lateral_shift(case_name: str):
    """ Define lateral shift for the plot
    
    This is used because spreading models may localize to the left or right
    of the middle point of the model domain.

    Shift factors were manually defined from initial plotting results to
    approximately center the middle of the plot with zone of localization.

    """
    shift_dict = {
        'case0': 0.0,
        }
    if case_name in shift_dict:
        return shift_dict[case_name]
    return 0.0


def get_variable_plot_parameters(local_plot_type: str) -> dict:
    """ Define variable plot parameters.

    These are the plot parameters that are specific to the type of plot being
    generated.

    """
    dx = get_lateral_shift(get_model_case_name_for_plotting())
    general_plot_params_dict = {
        LocalPlotTypes.GENERAL_ZOOM_OUT.name:
            {
                'figsize': (32.0, 8.0),
                'title_fontsize': 18,
                'axis_title_fontsize': 14,
                'axis_labels_fontsize': 14,
                'axis_ticks_fontsize': 14,
                'contour_label_fontsize': 10,
                'colorbar_labels_fontsize': 12,
                'colorbar_ticks_fontsize': 9,
                'text_box_font_size': 12,
                'axis_fraction_for_color_bar': 0.009,
                'axis_fraction_for_color_bar_gap': 0.025,
                'aspect_ratio': 1.0,
                'dimensions': (0.0 + dx, 500.0 + dx, 0.0, 160.0),
                'dim_close_up': (0.0 + dx, 500.0 + dx, 0.0, 160.0),
                'xyspacing': (50, 10.0),
                'xy_location_contour_legend': (505.0 + dx, -5.0),
                },
        LocalPlotTypes.ZOOM_IN_TO_HIGHRES.name:
            {
                'figsize': (32.0, 8.0),
                'title_fontsize': 18,
                'axis_title_fontsize': 14,
                'axis_labels_fontsize': 14,
                'axis_ticks_fontsize': 14,
                'contour_label_fontsize': 10,
                'colorbar_labels_fontsize': 12,
                'colorbar_ticks_fontsize': 9,
                'text_box_font_size': 12,
                'axis_fraction_for_color_bar': 0.009,
                'axis_fraction_for_color_bar_gap': 0.025,
                'aspect_ratio': 1.0,
                'dimensions': (150.0 + dx, 350.0 + dx, 0.0, 40.0),
                'dim_close_up': (150.0 + dx, 350.0 + dx, 0.0, 40.0),
                'xyspacing': (10, 5.0),
                'xy_location_contour_legend': (320.0 + dx, 4.0),
                },
        LocalPlotTypes.ZOOM_IN_2_HIGHRES_HF_GRAV.name:
            {
                'figsize': (20.0, 10.0),
                'title_fontsize': 12,
                'axis_title_fontsize': 10,
                'axis_labels_fontsize': 10,
                'axis_ticks_fontsize': 10,
                'contour_label_fontsize': 5,
                'colorbar_labels_fontsize': 5,
                'colorbar_ticks_fontsize': 5,
                'text_box_font_size': 8,
                'axis_fraction_for_color_bar': 0.008,
                'axis_fraction_for_color_bar_gap': 0.025,
                'aspect_ratio': 1.0,
                'dimensions': (150.0 + dx, 350.0 + dx, 0.0, 40.0),
                'dim_close_up': (150.0 + dx, 350.0 + dx, 0.0, 40.0),
                'xyspacing': (10, 5.0),
                'xy_location_contour_legend': (320.0 + dx, 5.0),
                }
        }
    return general_plot_params_dict[local_plot_type]


def get_plot_parameters(
        model_case_name: str,
        use_close_up: bool = False
) -> dict:
    """ Define 2D model plot parameters.

    Inputs
    ------
    model_case_name: str
        The case name of the model used to define the output directory. See
        paths.py for more information.

    """
    case_inputs = get_case_inputs(model_case_name)

    (
        nsteps
    ) = input_preprocessing.calculate_number_of_output_time_steps_for_velocity_stepping(
        case_inputs,
        istart = get_istart_iend()[0],
        iend = get_istart_iend()[1]
        )

    materials_input_dict, lib_path = get_materials_input_dict()
    output_dir = get_output_path(
        model_case_name=model_case_name, drive_number_id=2)
    plot_params = {
        'plot_output_path': os.path.join(output_dir, 'plots'),
        'material_library_file_path': lib_path,
        'materials_input_dict': materials_input_dict,
        'model_output_path': output_dir,
        'nsteps': nsteps,
        'nskip': 0,
        'istart': get_istart_iend()[0],
        'iplot_contour_labels': 1,
        'color_map': 'gnuplot2',
        'figure_dpi': 150.0,
        'extension': '.pdf',
        'length_units': 'km',
        'time_units': 'Myr',
        'velocity_units': 'cm/yr',
        'use_close_up': use_close_up,
        'number_format': '%6.1f',
        'rightside_up': True,
        }
    variable_params = get_variable_plot_parameters(define_local_plot_type())
    plot_params.update(variable_params)
    return plot_params


def get_model_plots_2d(model_case_name: str) -> ModelPlots2D:
    """ Get model plots 2D.
    """
    plot_parameters = get_plot_parameters(model_case_name)
    mp2d = ModelPlots2D(**plot_parameters)
    return mp2d


def plot_scalars(model_case_name: str) -> None:
    """ Plot scalars.
    """
    mp2d = get_model_plots_2d(model_case_name)
    scalar_names = ScalarNames()

    mp2d.grid_plots.plot_scalars(
        scalar_names.viscosity,
        contour_interval=0.5,
        minimum_value=18.0,
        maximum_value=26.0,
        plot_contours=False,
        grid_plot_type='nomesh'
        )
    use_all = False
    if use_all:
        mp2d.grid_plots.plot_scalars(
            scalar_names.shear_stress,
            contour_interval=5,
            minimum_value=-40,
            maximum_value=40,
            plot_contours=True,
            grid_plot_type='nomesh'
            )

        mp2d.grid_plots.plot_scalars(
            scalar_names.shear_plastic_failure,
            contour_interval=0.1,
            minimum_value=0,
            maximum_value=1,
            plot_contours=False,
            grid_plot_type='nomesh'
            )
        mp2d.grid_plots.plot_scalars(
            scalar_names.strainrate,
            contour_interval=0.5,
            minimum_value=-18.0,
            maximum_value=-13.0,
            plot_contours=False,
            grid_plot_type='nomesh'
            )
        mp2d.grid_plots.plot_scalars(
            scalar_names.thermal_conductivity,
            contour_interval=0.5,
            minimum_value=0,
            maximum_value=20,
            plot_contours=True,
            grid_plot_type='nomesh'
            )
        mp2d.grid_plots.plot_scalars(
            scalar_names.temperature,
            contour_interval=100.0,
            minimum_value=0.01,
            maximum_value=1400.0,
            #excluded_value=0.0,
            plot_contours=True,
            grid_plot_type='nomesh'
            )
        mp2d.grid_plots.plot_scalars(
            scalar_names.density,
            contour_interval=50,
            minimum_value=2700,
            maximum_value=3300,
            plot_contours=True,
            grid_plot_type='nomesh'
            )
        mp2d.grid_plots.plot_scalars(
            scalar_names.velocity_x,
            contour_interval=0.01,
            minimum_value=-0.1,
            maximum_value=0.1,
            plot_contours=False,
            grid_plot_type='nomesh'
            )
        mp2d.grid_plots.plot_scalars(
            scalar_names.velocity_y,
            contour_interval=0.01,
            minimum_value=-0.1,
            maximum_value=0.1,
            plot_contours=False,
            grid_plot_type='nomesh'
            )
        mp2d.grid_plots.plot_scalars(
            scalar_names.shear_stress,
            contour_interval=5,
            minimum_value=-40,
            maximum_value=40,
            plot_contours=True,
            grid_plot_type='nomesh'
            )
        mp2d.grid_plots.plot_scalars(
            scalar_names.shear_plastic_failure,
            contour_interval=0.1,
            minimum_value=0,
            maximum_value=1,
            plot_contours=False,
            grid_plot_type='nomesh'
            )
        mp2d.grid_plots.plot_scalars(
            scalar_names.velocity_mag,
            contour_interval=0.5,
            minimum_value=0,
            maximum_value=5.0,
            plot_contours=False,
            grid_plot_type='nomesh'
            )
        mp2d.grid_plots.plot_scalars(
            scalar_names.viscosity,
            contour_interval=0.5,
            minimum_value=18.0,
            maximum_value=26.0,
            plot_contours=True,
            grid_plot_type='nomesh'
            )


def plot_velocity(model_case_name: str) -> None:
    """ Plot velocity.
    """
    mp2d = get_model_plots_2d(model_case_name)
    mp2d.grid_plots.plot_velocity()


def plot_markers(model_case_name: str) -> None:
    """ Plot markers.
    """
    case_inputs = get_case_inputs(model_case_name)
    keys = CaseInputsKeys()

    eruption_internal_myr = case_inputs[keys.eruption_interval_yr]/1e6

    mp2d = get_model_plots_2d(model_case_name)

    plot_types = MarkerPlotNames
    cmap_names = MatplotlibColorMapNames

    local_plot_type = define_local_plot_type()

    mp2d.marker_plots.set_parameters(
        marker_size=get_marker_size(local_plot_type),
        plot_contour_labels=True,
        contour_line_color='black', decimation_factor_scatter_overlay=4,
        ny_contour_grid=50, nx_contour_grid=100, colorbar_shift_factor=4.0,
        decimation_factor=1
        )
    mp2d.marker_plots.set_parameters(
        plot_topography=True
        )
    mp2d.marker_plots.set_parameters(
        plot_base_level=False, base_level_line_color='blue',
        base_level_line_width=0.5,
        )
    mp2d.marker_plots.set_parameters(
        plot_meltfrac_contours=True, melt_fraction_min=0.0,
        melt_fraction_max=0.2, melt_fraction_contour_interval=0.05,
        meltfrac_contour_color='red', meltfrac_number_format='%6.2f',
        meltfrac_label_rightside_up=True
        )
    mp2d.marker_plots.set_parameters(
        plot_plastic_strain=True, strain_min=1.0, strain_max=6.0,
        strain_contour_interval=0.25, strain_cmap=cmap_names.inferno
        )
    mp2d.marker_plots.set_parameters(
        plot_plastic_strain_rate=True, strain_rate_min=-18, strain_rate_max=-12,
        strain_rate_contour_interval=0.5, strain_rate_cmap=cmap_names.Reds
        )
    mp2d.marker_plots.set_parameters(
        plot_sediment_age=False,
        age_min=19, age_max=26.0,
        age_contour_interval=0.25, age_cmap='hsv_r'
        )
    mp2d.marker_plots.set_parameters(
        plot_volcanics_age=True,
        age_min_volcanics=0, age_max_volcanics=15.0,
        age_contour_interval_volcanics=input_preprocessing.get_age_contour_interval_volcanics(
            eruption_internal_myr)
        )
    mp2d.marker_plots.set_parameters(
        plot_intrusive_age=True,
        age_min_intrusive=0, age_max_intrusive=15.0,
        age_contour_interval_intrusive=0.2,
        )
    mp2d.marker_plots.set_parameters(
        plot_serpentinization=False,
        serpentinization_min=0.0, serpentinization_max=1.0,
        serpentinization_contour_interval=0.1,
        serpentinization_cmap=cmap_names.Greens
        )
    mp2d.marker_plots.set_parameters(
        plot_temperature_contours=True,
        temperature_min=100.0, temperature_max=1500.0,
        temperature_contour_interval=100.0, temperature_contour_color='black'
        )
    mp2d.marker_plots.set_parameters(
        heatflow_min=25.0, heatflow_max=400, heatflow_spacing = 50,
        gravity_min=-300.0, gravity_max=300.0, gravity_spacing=50,
        subplot_spacing_fraction=0.3, height_ratios=[1, 1.25, 3]
        )

    if local_plot_type == LocalPlotTypes.ZOOM_IN_2_HIGHRES_HF_GRAV.name:
        mp2d.marker_plots.plot_markers(
            plot_type=plot_types.CompositionHeatFlowGravity.name)
    else:
        mp2d.marker_plots.plot_markers(
            plot_type=plot_types.Composition.name)


def plot_yield_stress(model_case_name: str) -> None:
    """ Plot yield stress.
    """
    mp2d = get_model_plots_2d(model_case_name)
    mp2d.rheology_plots.materials.override_material_properties(
        **get_case_inputs(model_case_name))
    mp2d.rheology_plots.plot_yield_strength(
        depth_plot_limit=150.0,
        depth_axis_spacing=10.0,
        temperature_plot_limit=1400.0,
        temperature_axis_spacing=200.0,
        yield_stress_plot_limit=600.0,
        yield_stress_axis_spacing=50.0,
        strain_rate=1e-15,
        thickness_upr_cont_crust_meters=22_000.0,
        thickness_lwr_cont_crust_meters=13_000.0,
        thickness_lithosphere_meters=125_000.0,
        thickness_asthenosphere_meters=25_000.0,
        dy_meters=100.0,
        expansivity=3e-5,
        compressibility=1e-11,
        density_upper_continental_crust=2822.72,
        density_lower_continental_crust=2900.0,
        density_mantle_lithosphere=3250.0,
        density_asthenosphere=3300.0,
        use_linear_segments=False,
        temperature_top_celsius=0.0,
        temperature_moho_celsius=600.0,
        temperature_base_lith_celsius=1330.0,
        adiabatic_gradient_kelvin_km=0.4,
        conductivity_upper_crust=2.25,
        conductivity_lower_crust=2.2,
        conductivity_mantle=2.0,
        heat_production_upper_crust=1.8e-6,
        heat_production_lower_crust=0.5e-6,
        heat_production_mantle=0.0,
        thickness_thermal_lithosphere=125_000.0,
        iuse_fluid_pressure_for_yield=1,
        extension='.pdf',
        )


def plot_stokes_convergence(model_case_name: str):
    """ Plot stokes convergence.
    """
    mp2d = get_model_plots_2d(model_case_name)
    mp2d.make_stokes_convergence_plot(
        figsize=(16, 4),
        xmin=0,
        xmax=10000,
        xspacing=500,
        log_l2_ymin=-7,
        log_l2_ymax=1,
        log_l2_yspacing=1,
        plot_title=model_case_name,
        axis_label_size=8,
        axis_title_size=8,
        legend_font_size=8,
        xtick_label_size=6,
        ytick_label_size=6,
        annotation_font_size=4
        )


def _plot_magmatic_crust_height(model_case_name: str):
    output_directory_path = get_output_path(model_case_name=model_case_name)
    plot_magmatic_crust_height(
        output_directory_path,
        model_case_name,
        xdim=(0, 20),
        ydim=(0, 10000),
        yticks=(0, 10000, 500)
        )


def get_func(option_name: str) -> Callable[[None], None]:
    """ Get function based on option name.
    """
    func_dict = {
        'plot_markers': plot_markers,
        'plot_scalars': plot_scalars,
        'plot_velocity': plot_velocity,
        'plot_yield_stress': plot_yield_stress,
        'plot_stokes_convergence': plot_stokes_convergence,
        'plot_magmatic_crust_height': _plot_magmatic_crust_height
        }
    if option_name not in func_dict:
        raise ValueError(
            f'{option_name} is not a valid option. '
            f'Valid option are {list(func_dict.keys())}'
            )
    return func_dict[option_name]


if __name__ == '__main__':
    get_func(sys.argv[1])(model_case_name=get_model_case_name_for_plotting())
