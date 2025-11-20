""" Plotting functions.

Usage
-----
To plot markers:
    % python plot.py plot_markers
To plot scalars:
    % python plot.py plot_scalars
To plot velocity:
    % python plot.py plot_velocity

"""
import os
import sys
from typing import Callable
from materials import get_materials_input_dict
from paths import get_output_path
from earthbox.plot_tools import ModelPlots2D
from earthbox.plot_tools import ScalarNames
from earthbox.plot_tools import MarkerPlotNames
from earthbox.plot_tools.utils.color_maps import MatplotlibColorMapNames


def get_plot_parameters(use_close_up=False) -> dict:
    """ Define 2D model plot parameters.
    """
    materials_input_dict, lib_path = get_materials_input_dict()
    output_dir = get_output_path()
    return {
        'plot_output_path': os.path.join(output_dir, 'plots'),
        'material_library_file_path': lib_path,
        'materials_input_dict': materials_input_dict,
        'model_output_path': output_dir,
        'nsteps': 200,
        'nskip': 0,
        'istart': 0,
        'iplot_contour_labels': 1,
        'color_map': 'bwr',
        'figure_dpi': 150.0,
        'figsize': (20.0, 10.0),
        'aspect_ratio': get_aspect_ration(use_close_up),
        'length_units': 'km',
        'time_units': 'Myr',
        'velocity_units': 'cm/yr',
        'axis_fraction_for_color_bar': 0.008,
        'axis_fraction_for_color_bar_gap': 0.008,
        'dimensions': (0.0, 1000.0, 0.0, 160.0),
        'use_close_up': use_close_up,
        'dim_close_up': (200.0, 600.0, 0.0, 50.0),
        'xyspacing': get_xyspacing(use_close_up),
        'xy_location_contour_legend': get_xy_location_contour_legend(use_close_up),
        'title_fontsize': 12,
        'axis_title_fontsize': 10,
        'axis_labels_fontsize': 10,
        'axis_ticks_fontsize': 10,
        'contour_label_fontsize': 5,
        'colorbar_ticks_fontsize': 4
        }


def get_xyspacing(use_close_up: bool) -> tuple:
    """ Get x and y spacing.
    """
    xy_spacing = (50.0, 10.0)
    if use_close_up is True:
        xy_spacing = (25.0, 5.0)
    return xy_spacing


def get_xy_location_contour_legend(use_close_up: bool) -> tuple:
    """ Get xy location contour legend.
    """
    xy_location = (705.0, -2.5)
    if use_close_up is True:
        xy_location = (1005.0, -5.0)
    return xy_location


def get_aspect_ration(use_close_up: bool) -> float:
    """ Get aspect ratio.
    """
    aspect_ratio = 1.5
    if use_close_up is True:
        aspect_ratio = 4.0
    return aspect_ratio


def get_model_plots_2d() -> ModelPlots2D:
    """ Get model plots 2D.
    """
    plot_parameters = get_plot_parameters()
    mp2d = ModelPlots2D(**plot_parameters)
    return mp2d


def plot_scalars():
    """ Plot scalars.
    """
    mp2d = get_model_plots_2d()
    scalar_names = ScalarNames()
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
        mp2d.grid_plots.plot_scalars(
            scalar_names.strainrate,
            contour_interval=0.5,
            minimum_value=-18.0,
            maximum_value=-13.0,
            plot_contours=False,
            grid_plot_type='nomesh'
            )


def plot_velocity():
    """ Plot velocity.
    """
    mp2d = get_model_plots_2d()
    mp2d.grid_plots.plot_velocity()


def plot_markers():
    """ Plot markers.
    """
    mp2d = get_model_plots_2d()
    plot_types = MarkerPlotNames
    cmap_names = MatplotlibColorMapNames
    mp2d.marker_plots.set_parameters(
        marker_size=0.3, plot_contour_labels=True,
        contour_line_color='black', decimation_factor_scatter_overlay=4,
        ny_contour_grid=50, nx_contour_grid=100, colorbar_shift_factor=4.0
        )
    mp2d.marker_plots.set_parameters(
        plot_topography=True
        )
    mp2d.marker_plots.set_parameters(
        plot_base_level=True, base_level_line_color='blue',
        base_level_line_width=0.5,
        )
    mp2d.marker_plots.set_parameters(
        plot_meltfrac_contours=True, melt_fraction_min=0.0,
        melt_fraction_max=0.2, melt_fraction_contour_interval=0.1,
        meltfrac_contour_color='red'
        )
    mp2d.marker_plots.set_parameters(
        plot_plastic_strain=True, strain_min=1.0, strain_max=6.0,
        strain_contour_interval=0.25, strain_cmap=cmap_names.inferno
        )
    mp2d.marker_plots.set_parameters(
        plot_sediment_age=True, age_min=0.0, age_max=90.0,
        age_contour_interval=2.5, age_cmap='hsv_r'
        )
    mp2d.marker_plots.set_parameters(
        plot_serpentinization=True,
        serpentinization_min=0.0, serpentinization_max=1.0,
        serpentinization_contour_interval=0.1,
        serpentinization_cmap=cmap_names.Greens
        )
    mp2d.marker_plots.set_parameters(
        plot_temperature_contours=True,
        temperature_min=0.01, temperature_max=1400.0,
        temperature_contour_interval=100.0, temperature_contour_color='white'
        )
    mp2d.marker_plots.set_parameters(
        heatflow_min=25.0, heatflow_max=200, heatflow_spacing = 25,
        gravity_min=-300.0, gravity_max=300.0, gravity_spacing=50,
        subplot_spacing_fraction=0.3, height_ratios=[1, 1.25, 3]
        )
    mp2d.marker_plots.plot_markers(plot_type=plot_types.CompositionHeatFlowGravity.name)


def plot_yield_stress():
    """ Plot yield stress.
    """
    mp2d = get_model_plots_2d()
    mp2d.rheology_plots.plot_yield_strength(
        depth_plot_limit=140.0,
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
        )


def plot_stokes_convergence():
    """ Plot stokes convergence.
    """
    mp2d = get_model_plots_2d()
    mp2d.make_stokes_convergence_plot(
        figsize=(80, 10),
        xspacing=5,
        log_l2_ymin=-5,
        log_l2_ymax=1,
        log_l2_yspacing=1,
        plot_title='Case0.1',
        axis_label_size=20,
        axis_title_size=20,
        legend_font_size=15,
        xtick_label_size=14,
        ytick_label_size=14,
        annotation_font_size=4
        )


def get_func(option_name: str) -> Callable[[None], None]:
    """ Get function based on option name.
    """
    func_dict = {
        'plot_markers': plot_markers,
        'plot_scalars': plot_scalars,
        'plot_velocity': plot_velocity,
        'plot_yield_stress': plot_yield_stress,
        'plot_stokes_convergence': plot_stokes_convergence
        }
    if option_name not in func_dict:
        raise ValueError(
            f'{option_name} is not a valid option. '
            f'Valid option are {list(func_dict.keys())}'
            )
    return func_dict[option_name]


if __name__ == '__main__':
    get_func(sys.argv[1])()
