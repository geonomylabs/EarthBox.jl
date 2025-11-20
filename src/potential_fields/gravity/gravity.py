#!/usr/bin/env python3
# -*- coding: utf-8 -*-
""" Functions used calculate gravity anomaly
"""
import os
import math
import numpy as np
import numpy.typing as npt
import numba
from numba import njit
from earthbox.model import PyModelData
from earthbox.utilities.printfuncs import timeit_memit
from earthbox.utilities.conversion_funcs \
    import meters_per_second_squared_to_milligal
from earthbox.grids.basic.spacing import get_basic_grid_spacing_arrays
from .cell.cell import calculate_gravity_anomaly_of_cell


@timeit_memit('Finished calculating gravity anomaly')
def calculate_gravity_anomaly(
        pymodel: PyModelData, # type: ignore
        output_dir_path: str
) -> None:
    """ Calculate Free Air and Bouguer gravity anomalies.
    """
    iuse_topo = pymodel.topography.parameters.model_options.iuse_topo.value
    if iuse_topo == 1:
        (
            topo_gridx, gravity_grid_mgal, gravity_grid_free_air_mgal
        ) = calculate_free_air_and_bouguer_active_model(pymodel.numba_model)

        save_gravity_grids_to_file(
            topo_gridx, gravity_grid_mgal, gravity_grid_free_air_mgal,
            output_dir_path
            )
    else:
        print(
            '!!! WARNING !!! Gravity can only be calculated if topography is '
            ' activated.'
            )


@njit
def calculate_free_air_and_bouguer_active_model(
        model: PyModelData # type: ignore
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """ Calculate gravity anomalies.

    This function calculates the gravity anomalies for Free Air and Bouguer
    corrections.

    Returns
    -------
    topo_gridx : npt.NDArray[np.float64]
        x-coordinates of topography grid nodes (meters).

    gravity_grid_mgal : npt.NDArray[np.float64]
        Gravity anomaly with Free Air correction applied (mgal).

    gravity_grid_bouguer_mgal : npt.NDArray[np.float64]
        Gravity anomaly with Bouguer correction applied (mgal).

    """
    gravity_constant = 6.6732e-11
    crustal_density = 2800.0 # kg/m^3, used for Bouguer correction

    bgridx = model.grids.arrays.basic.bgridx.array
    bgridy = model.grids.arrays.basic.bgridy.array
    xstp = model.grids.arrays.basic.xstp1.array
    ystp = model.grids.arrays.basic.ystp1.array

    rho_grid = model.stokes_continuity.arrays.density.rho1.array

    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    base_level_shift = model.topography.parameters.sealevel.base_level_shift.value
    # Use global base level
    y_sealevel = y_sealevel - base_level_shift

    gridt = model.topography.arrays.gridt.array

    topo_gridx = gridt[0, :].copy()
    gravity_grid_mgal = gravity_anomaly_loop_left_edge(
        y_sealevel, bgridx, bgridy, xstp, ystp, rho_grid,
        topo_gridx, gravity_constant
        )

    topo_gridy = gridt[1, :].copy()
    gravity_grid_free_air_mgal = calculate_free_air_gravity(
        y_sealevel, gravity_constant, crustal_density,
        topo_gridx, topo_gridy, gravity_grid_mgal
        )
    return topo_gridx, gravity_grid_mgal, gravity_grid_free_air_mgal


@timeit_memit('Finished calculating gravity anomaly')
@njit
def calculate_free_air_and_bouguer(
        bgridx: npt.NDArray[np.float64],
        bgridy: npt.NDArray[np.float64],
        rho_grid: npt.NDArray[np.float64],
        y_sealevel: float,
        topo_gridx: npt.NDArray[np.float64],
        topo_gridy: npt.NDArray[np.float64],
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """ Calculate gravity anomalies.

    This function calculates the gravity anomalies for Free Air and Bouguer
    corrections.

    Inputs
    ------
    bgridx : Array((xnum), dtype=np.float64)
        x-coordinates of grid nodes (meters).

    bgridy : Array((ynum), dtype=np.float64)
        y-coordinates of grid nodes (meters).

    rho_grid : Array((ynum, xnum), dtype=np.float64)
        Density of grid nodes (kg/m^3).

    y_sealevel : float
        y-coordinate of sealevel (meters).

    topo_gridx : Array((ntopo), dtype=np.float64)
        x-coordinates of topography grid nodes (meters).

    topo_gridy : Array((ntopo), dtype=np.float64)
        y-coordinates of topography grid nodes (meters).

    Returns
    -------
    gravity_grid_mgal : Array((ntopo), dtype=np.float64)
        Gravity anomaly with Free Air correction applied (mgal).

    gravity_grid_bouguer_mgal : Array((ntopo), dtype=np.float64)
        Gravity anomaly with Bouguer correction applied (mgal).

    """
    gravity_constant = 6.6732e-11
    crustal_density = 2800.0 # kg/m^3, used for Bouguer correction

    xstp, ystp = get_basic_grid_spacing_arrays(bgridx, bgridy)

    gravity_grid_mgal = gravity_anomaly_loop_left_edge(
        y_sealevel, bgridx, bgridy, xstp, ystp, rho_grid,
        topo_gridx, gravity_constant
        )

    gravity_grid_free_air_mgal = calculate_free_air_gravity(
        y_sealevel, gravity_constant, crustal_density,
        topo_gridx, topo_gridy, gravity_grid_mgal
        )
    return gravity_grid_mgal, gravity_grid_free_air_mgal


@njit(parallel=True)
def gravity_anomaly_loop_left_edge(
        y_sealevel: float,
        bgridx: npt.NDArray[np.float64],
        bgridy: npt.NDArray[np.float64],
        xstp: npt.NDArray[np.float64],
        ystp: npt.NDArray[np.float64],
        rho_grid: npt.NDArray[np.float64],
        topo_gridx: npt.NDArray[np.float64],
        gravity_constant: float
) -> npt.NDArray[np.float64]:
    """ Loop over grid to calculate gravity anomaly.

    This loop calculate free air gravity anomaly offshore and Bouguer anomaly
    onshore relative to sealevel.

    """
    # pylint: disable=too-many-arguments
    # pylint: disable=too-many-locals
    toponum = topo_gridx.shape[0]
    gravity_grid_mgal = np.zeros(toponum, dtype=np.float64)

    xnum = bgridx.shape[0]
    ynum = bgridy.shape[0]

    for itopo in numba.prange(toponum): # pylint: disable=not-an-iterable
        x_observer = topo_gridx[itopo]
        for i in range(ynum):
            y_grid = bgridy[i]
            # Only consider mass anomalies below sealevel
            if y_grid >= y_sealevel:
                for j in range(xnum):
                    x_grid = bgridx[j]
                    (
                        delta_density
                    ) = calculate_delta_density_relative_to_left_edge(
                        rho_grid, i, j)

                    (
                        x_upper_left_cell, cell_width
                    ) = calculate_horizontal_cell_geometry(
                        j, xstp, x_grid)

                    (
                        y_upper_left_cell, cell_height
                    ) = calculate_vertical_cell_geometry(
                        i, ystp, y_grid, y_sealevel)

                    if y_upper_left_cell >= y_sealevel:

                        (
                            gravity_anomaly
                        ) = calculate_gravity_anomaly_of_cell(
                            gravity_constant, delta_density,
                            x_upper_left_cell, y_upper_left_cell,
                            x_observer, y_sealevel, cell_height, cell_width,
                            beta=0.0
                            )

                        (
                            gravity_anomaly_mgal
                        ) = meters_per_second_squared_to_milligal(
                            gravity_anomaly)

                        gravity_grid_mgal[itopo] = (
                            gravity_grid_mgal[itopo]
                            + gravity_anomaly_mgal
                            )

    return gravity_grid_mgal


@njit
def calculate_free_air_gravity(
        y_sealevel: float,
        gravity_constant: float,
        crustal_density: float,
        topo_gridx: npt.NDArray[np.float64],
        topo_gridy: npt.NDArray[np.float64],
        gravity_grid_mgal: npt.NDArray[np.float64]
) -> npt.NDArray[np.float64]:
    """ Apply Bouguer correction.

    Returns
    -------
    gravity_grid_free_air_mgal : npt.NDArray[np.float64]
        Gravity anomaly with Bouguer correction applied yielding the Free Air
        gravity anomaly (mgal).

    """
    # pylint: disable=too-many-arguments
    nnodes = gravity_grid_mgal.shape[0]
    gravity_grid_free_air_mgal = np.zeros(nnodes, dtype=np.float64)

    y_topo_left_edge = topo_gridy[0]
    topo_height_left_edge = max(0.0, y_sealevel - y_topo_left_edge)

    toponum = topo_gridx.shape[0]
    for itopo in range(toponum):
        y_topo = topo_gridy[itopo]
        if y_topo < y_sealevel:
            topo_height = max(0.0, y_sealevel - y_topo)
            delta_topo_height = max(0.0, topo_height - topo_height_left_edge)
            correction = (
                2.0*math.pi*gravity_constant*crustal_density*delta_topo_height)
            correction_mgal = meters_per_second_squared_to_milligal(correction)
            gravity_grid_free_air_mgal[itopo] = (
                gravity_grid_mgal[itopo] + correction_mgal)
        else:
            gravity_grid_free_air_mgal[itopo] = gravity_grid_mgal[itopo]
    return gravity_grid_free_air_mgal


@njit
def calculate_delta_density_relative_to_left_edge(
        rho_grid: npt.NDArray[np.float64],
        i: int,
        j: int
) -> float:
    """ Get delta density relative to left edge of grid.
    """
    delta_density = rho_grid[i, j] - rho_grid[i, 0]
    return delta_density


def save_gravity_grids_to_file(
        topo_gridx: npt.NDArray[np.float64],
        gravity_grid_mgal: npt.NDArray[np.float64],
        gravity_grid_free_air_mgal: npt.NDArray[np.float64],
        output_dir_path: str
) -> None:
    """ Save gravity grids to file.
    """
    filename = os.path.join(output_dir_path, 'gravity_grids.txt')
    with open(filename, 'w', encoding='utf-8') as file:
        file.write('Topography_x_(m) Bouguer_(mgal) Free_Air_(mgal)\n')
        for i, x in enumerate(topo_gridx):
            gravity_bouguer = gravity_grid_mgal[i]
            gravity_free_air = gravity_grid_free_air_mgal[i]
            file.write(f'{x} {gravity_bouguer} {gravity_free_air}\n')


@njit
def calculate_vertical_cell_geometry(
        i: int,
        ystp: npt.NDArray[np.float64],
        y_grid: float,
        y_sealevel: float
) -> tuple[float, float]:
    """ Calculate vertical cell geometry taking sealevel into account.

    Inputs
    ------
    i : int
        Index of grid node in y-direction.

    ystp : npt.NDArray[np.float64]
        y-spacing of grid nodes (meters).

    y_grid : float
        y-coordinate of grid node (meters).

    y_sealevel : float
        y-coordinate of sealevel (meters).

    Returns
    -------
    y_upper_left_cell : float
        y-coordinate of upper left corner of cell (meters).

    cell_height : float
        Height of cell (meters).

    """
    cell_height = get_cell_size_y(ystp, i)
    y_upper_left_cell, truncated_cell_height = get_y_upper_left(
        y_grid, cell_height, y_sealevel)
    # Correct cell height if truncated
    cell_height = cell_height - truncated_cell_height
    return y_upper_left_cell, cell_height


@njit
def get_x_upper_left(x_grid: float, cell_width) -> float:
    """ Get x-coordinate of upper left corner of cell.
    """
    return x_grid - cell_width/2.0


@njit
def get_y_upper_left(y_grid: float, cell_height, y_sealevel: float) -> float:
    """ Get y-coordinate of upper left corner of cell.

    The cell is truncated by the sealevel.

    """
    y_upper_left = y_grid - cell_height/2.0

    if y_upper_left < y_sealevel:
        truncated_cell_height = y_sealevel - y_upper_left
    else:
        truncated_cell_height = 0.0

    y_upper_left = max(y_upper_left, y_sealevel)
    return y_upper_left, truncated_cell_height


@njit
def calculate_horizontal_cell_geometry(
        j: int,
        xstp: npt.NDArray[np.float64],
        x_grid: float
) -> tuple[float, float]:
    """ Calculate horizontal cell geometry.

    Inputs
    ------
    j : int
        Index of grid node in x-direction.

    xstp : npt.NDArray[np.float64]
        x-spacing of grid nodes (meters).

    x_grid : float
        x-coordinate of grid node (meters).

    Returns
    -------
    x_upper_left_cell : float
        x-coordinate of upper left corner of cell (meters).

    cell_width : float
        Width of cell (meters).

    """
    cell_width = get_cell_size_x(xstp, j)
    x_upper_left_cell = get_x_upper_left(x_grid, cell_width)
    return x_upper_left_cell, cell_width


@njit
def get_cell_size_x(
        xstp: npt.NDArray[np.float64],
        j: int
) -> float:
    """ Get cell size in x-direction.
    """
    nstep = xstp.shape[0]
    if j == 0:
        cell_size_x = xstp[j]
    elif j < nstep - 1:
        cell_size_x = (xstp[j-1] + xstp[j])/2.0
    else:
        cell_size_x = xstp[nstep-1]
    return cell_size_x


@njit
def get_cell_size_y(
        ystp: npt.NDArray[np.float64],
        i: int
) -> float:
    """ Get cell size in y-direction.
    """
    nstep = ystp.shape[0]
    if i == 0:
        cell_size_y = ystp[i]
    elif i < nstep - 1:
        cell_size_y = (ystp[i-1] + ystp[i])/2.0
    else:
        cell_size_y = ystp[nstep-1]
    return cell_size_y
