#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Script for working with models where parameters are read from input files.

Usage
-----
To run model:
    % python read_run.py run_model
To plot markers:
    % python read_run.py plot_markers
To plot scalars:
    % python read_run.py plot_scalars
To plot velocity:
    % python read_run.py plot_velocity

Note that plots can be made in time loops by activating loop plotting using
the setup_loop_plotting function.

Input Files
-----------
Model input file: model.yml
    File contains lists with value, units and description for each model
    parameter associated with the model.

Materials input file: materials.yml
    This file contains the integer ID of each material followed by material
    name from library, material type, material domain and rgb colors.

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


def run_model():
    """ Run model.
    """
    eb = EarthBox(
        paths={
            'model_input_file': "model.yml",
            'materials_input_file': "materials.yml",
            'materials_library_file': get_materials_lib_path(),
            'output_dir': get_output_path()
            },
        )
    eb.model_manager.config.solver.use_mumps = False
    eb.model_manager.config.solver.nproc = 1
    eb.model_manager.initialize_model()
    eb.run_time_steps(make_backup=False)


def get_model_plots_2d():
    """ Get model plots 2D.
    """
    mp2d = ModelPlots2D(
        plot_output_path=os.path.join(get_output_path(), 'plots'),
        material_library_file_path=get_materials_lib_path(),
        material_model_file_path='materials.yml',
        model_output_path=get_output_path(),
        nsteps=100,
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
        xyspacing=(2.0, 1.0)
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
        marker_size=0.25
        )


def get_func(option_name: str) -> Callable[[None], None]:
    """ Get function based on option name.
    """
    func_dict = {
        'run_model': run_model,
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
