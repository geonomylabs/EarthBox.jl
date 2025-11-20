module Gravity

include("gravity_cell/GravityCell.jl")

import EarthBox.ModelDataContainer: ModelData
import .GravityCell: calculate_gravity_anomaly_of_cell
import EarthBox.StaggeredGrid.Spacing: get_basic_grid_spacing_arrays!
import EarthBox.ConversionFuncs: meters_per_second_squared_to_milligal

""" Calculate Free Air and Bouguer gravity anomalies.

# Arguments
- `model::ModelData`: The model data container
- `output_dir_path::String`: Path to output directory

# Returns
- `Nothing`
"""
function calculate_gravity_anomaly(
    model::ModelData,
    output_dir_path::String
)::Nothing
    iuse_topo = model.topography.parameters.model_options.iuse_topo.value
    if iuse_topo == 1
        (
            topo_gridx, gravity_grid_mgal, gravity_grid_free_air_mgal
        ) = calculate_free_air_and_bouguer_active_model(model)

        save_gravity_grids_to_file(
            topo_gridx, gravity_grid_mgal, gravity_grid_free_air_mgal,
            output_dir_path
        )
    else
        println(
            "!!! WARNING !!! Gravity can only be calculated if topography is "
            * "activated."
        )
    end
    return nothing
end

""" Calculate gravity anomalies.

This function calculates the gravity anomalies for Free Air and Bouguer
corrections.

# Arguments
- `model::ModelData`: The model data container

# Returns
- `topo_gridx::Vector{Float64}`: x-coordinates of topography grid nodes (meters)
- `gravity_grid_mgal::Vector{Float64}`: Gravity anomaly with Free Air correction 
  applied (mgal)
- `gravity_grid_free_air_mgal::Vector{Float64}`: Gravity anomaly with Bouguer 
  correction applied (mgal)
"""
function calculate_free_air_and_bouguer_active_model(
    model::ModelData
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    gravity_constant = 6.6732e-11
    crustal_density = 2800.0  # kg/m^3, used for Bouguer correction

    bgridx = model.grids.arrays.basic.gridx_b.array
    bgridy = model.grids.arrays.basic.gridy_b.array
    xstp, ystp = get_basic_grid_spacing_arrays!(bgridx, bgridy)

    rho_grid = model.stokes_continuity.arrays.density.rho1.array

    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    base_level_shift = model.topography.parameters.sealevel.base_level_shift.value
    # Use global base level
    y_sealevel = y_sealevel - base_level_shift

    gridt = model.topography.arrays.gridt.array

    topo_gridx = gridt[1, :]
    gravity_grid_mgal = gravity_anomaly_loop_left_edge(
        y_sealevel, bgridx, bgridy, xstp, ystp, rho_grid,
        topo_gridx, gravity_constant
    )

    topo_gridy = gridt[2, :]
    gravity_grid_free_air_mgal = calculate_free_air_gravity(
        y_sealevel, gravity_constant, crustal_density,
        topo_gridx, topo_gridy, gravity_grid_mgal
    )
    return topo_gridx, gravity_grid_mgal, gravity_grid_free_air_mgal
end

""" Calculate gravity anomalies.

This function calculates the gravity anomalies for Free Air and Bouguer
corrections.

# Arguments
- `bgridx::Vector{Float64}`: x-coordinates of grid nodes (meters)
- `bgridy::Vector{Float64}`: y-coordinates of grid nodes (meters)
- `rho_grid::Matrix{Float64}`: Density of grid nodes (kg/m^3)
- `y_sealevel::Float64`: y-coordinate of sealevel (meters)
- `topo_gridx::Vector{Float64}`: x-coordinates of topography grid nodes (meters)
- `topo_gridy::Vector{Float64}`: y-coordinates of topography grid nodes (meters)

# Returns
- `gravity_grid_mgal::Vector{Float64}`: Gravity anomaly with Free Air correction 
  applied (mgal)
- `gravity_grid_free_air_mgal::Vector{Float64}`: Gravity anomaly with Bouguer 
  correction applied (mgal)
"""
function calculate_free_air_and_bouguer(
    bgridx::Vector{Float64},
    bgridy::Vector{Float64},
    rho_grid::Matrix{Float64},
    y_sealevel::Float64,
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}}
    gravity_constant = 6.6732e-11
    crustal_density = 2800.0  # kg/m^3, used for Bouguer correction

    xstp, ystp = get_basic_grid_spacing_arrays!(bgridx, bgridy)

    t1 = time()
    gravity_grid_mgal = gravity_anomaly_loop_left_edge(
        y_sealevel, bgridx, bgridy, xstp, ystp, rho_grid,
        topo_gridx, gravity_constant
    )
    t2 = time()
    println("Time taken to calculate gravity anomaly: $(t2 - t1) seconds")

    t1 = time()
    gravity_grid_free_air_mgal = calculate_free_air_gravity(
        y_sealevel, gravity_constant, crustal_density,
        topo_gridx, topo_gridy, gravity_grid_mgal
    )
    t2 = time()
    println("Time taken to calculate free air gravity anomaly: $(t2 - t1) seconds")
    return gravity_grid_mgal, gravity_grid_free_air_mgal
end

""" Loop over grid to calculate gravity anomaly.

This loop calculate free air gravity anomaly offshore and Bouguer anomaly
onshore relative to sealevel.

# Arguments
- `y_sealevel::Float64`: y-coordinate of sealevel (meters)
- `bgridx::Vector{Float64}`: x-coordinates of grid nodes (meters)
- `bgridy::Vector{Float64}`: y-coordinates of grid nodes (meters)
- `xstp::Vector{Float64}`: x-spacing of grid nodes (meters)
- `ystp::Vector{Float64}`: y-spacing of grid nodes (meters)
- `rho_grid::Matrix{Float64}`: Density of grid nodes (kg/m^3)
- `topo_gridx::Vector{Float64}`: x-coordinates of topography grid nodes (meters)
- `gravity_constant::Float64`: Gravitational constant (m³/kg/s²)

# Returns
- `gravity_grid_mgal::Vector{Float64}`: Gravity anomaly (mgal)
"""
function gravity_anomaly_loop_left_edge(
    y_sealevel::Float64,
    bgridx::Vector{Float64},
    bgridy::Vector{Float64},
    xstp::Vector{Float64},
    ystp::Vector{Float64},
    rho_grid::Matrix{Float64},
    topo_gridx::Vector{Float64},
    gravity_constant::Float64
)::Vector{Float64}
    toponum = length(topo_gridx)
    gravity_grid_mgal = zeros(Float64, toponum)

    xnum = length(bgridx)
    ynum = length(bgridy)

    for itopo in 1:toponum
        x_observer = topo_gridx[itopo]
        for j in 1:xnum  # Swapped loop order for column-major arrays
            x_grid = bgridx[j]
            for i in 1:ynum
                y_grid = bgridy[i]
                # Only consider mass anomalies below sealevel
                if y_grid >= y_sealevel
                    delta_density = calculate_delta_density_relative_to_left_edge(
                        rho_grid, i, j)

                    (x_upper_left_cell, cell_width) = calculate_horizontal_cell_geometry(
                        j, xstp, x_grid)

                    (y_upper_left_cell, cell_height) = calculate_vertical_cell_geometry(
                        i, ystp, y_grid, y_sealevel)

                    if y_upper_left_cell >= y_sealevel
                        gravity_anomaly = calculate_gravity_anomaly_of_cell(
                            gravity_constant, delta_density,
                            x_upper_left_cell, y_upper_left_cell,
                            x_observer, y_sealevel, cell_height, cell_width,
                            beta=0.0
                        )

                        gravity_anomaly_mgal = meters_per_second_squared_to_milligal(
                            gravity_anomaly)

                        gravity_grid_mgal[itopo] = (
                            gravity_grid_mgal[itopo] + gravity_anomaly_mgal
                        )
                    end
                end
            end
        end
    end

    return gravity_grid_mgal
end

""" Apply Bouguer correction.

# Arguments
- `y_sealevel::Float64`: y-coordinate of sealevel (meters)
- `gravity_constant::Float64`: Gravitational constant (m³/kg/s²)
- `crustal_density::Float64`: Crustal density (kg/m^3)
- `topo_gridx::Vector{Float64}`: x-coordinates of topography grid nodes (meters)
- `topo_gridy::Vector{Float64}`: y-coordinates of topography grid nodes (meters)
- `gravity_grid_mgal::Vector{Float64}`: Gravity anomaly (mgal)

# Returns
- `gravity_grid_free_air_mgal::Vector{Float64}`: Gravity anomaly with Bouguer 
  correction applied yielding the Free Air gravity anomaly (mgal)
"""
function calculate_free_air_gravity(
    y_sealevel::Float64,
    gravity_constant::Float64,
    crustal_density::Float64,
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    gravity_grid_mgal::Vector{Float64}
)::Vector{Float64}
    nnodes = length(gravity_grid_mgal)
    gravity_grid_free_air_mgal = zeros(Float64, nnodes)

    y_topo_left_edge = topo_gridy[1]
    topo_height_left_edge = max(0.0, y_sealevel - y_topo_left_edge)

    toponum = length(topo_gridx)
    for itopo in 1:toponum
        y_topo = topo_gridy[itopo]
        if y_topo < y_sealevel
            topo_height = max(0.0, y_sealevel - y_topo)
            delta_topo_height = max(0.0, topo_height - topo_height_left_edge)
            correction = (
                2.0 * π * gravity_constant * crustal_density * delta_topo_height)
            correction_mgal = meters_per_second_squared_to_milligal(correction)
            gravity_grid_free_air_mgal[itopo] = (
                gravity_grid_mgal[itopo] + correction_mgal)
        else
            gravity_grid_free_air_mgal[itopo] = gravity_grid_mgal[itopo]
        end
    end
    return gravity_grid_free_air_mgal
end

""" Get delta density relative to left edge of grid.

# Arguments
- `rho_grid::Matrix{Float64}`: Density grid
- `i::Int`: y-index
- `j::Int`: x-index

# Returns
- `delta_density::Float64`: Density difference relative to left edge
"""
function calculate_delta_density_relative_to_left_edge(
    rho_grid::Matrix{Float64},
    i::Int,
    j::Int
)::Float64
    delta_density = rho_grid[i, j] - rho_grid[i, 1]
    return delta_density
end

""" Save gravity grids to file.

# Arguments
- `topo_gridx::Vector{Float64}`: x-coordinates of topography grid nodes (meters)
- `gravity_grid_mgal::Vector{Float64}`: Gravity anomaly (mgal)
- `gravity_grid_free_air_mgal::Vector{Float64}`: Free air gravity anomaly (mgal)
- `output_dir_path::String`: Path to output directory

# Returns
- `Nothing`
"""
function save_gravity_grids_to_file(
    topo_gridx::Vector{Float64},
    gravity_grid_mgal::Vector{Float64},
    gravity_grid_free_air_mgal::Vector{Float64},
    output_dir_path::String
)::Nothing
    filename = joinpath(output_dir_path, "gravity_grids.txt")
    open(filename, "w") do file
        write(file, "Topography_x_(m) Bouguer_(mgal) Free_Air_(mgal)\n")
        for i in 1:length(topo_gridx)
            x = topo_gridx[i]
            gravity_bouguer = gravity_grid_mgal[i]
            gravity_free_air = gravity_grid_free_air_mgal[i]
            write(file, "$x $gravity_bouguer $gravity_free_air\n")
        end
    end
    return nothing
end

""" Calculate vertical cell geometry taking sealevel into account.

# Arguments
- `i::Int`: Index of grid node in y-direction
- `ystp::Vector{Float64}`: y-spacing of grid nodes (meters)
- `y_grid::Float64`: y-coordinate of grid node (meters)
- `y_sealevel::Float64`: y-coordinate of sealevel (meters)

# Returns
- `y_upper_left_cell::Float64`: y-coordinate of upper left corner of cell (meters)
- `cell_height::Float64`: Height of cell (meters)
"""
function calculate_vertical_cell_geometry(
    i::Int,
    ystp::Vector{Float64},
    y_grid::Float64,
    y_sealevel::Float64
)::Tuple{Float64, Float64}
    cell_height = get_cell_size_y(ystp, i)
    y_upper_left_cell, truncated_cell_height = get_y_upper_left(
        y_grid, cell_height, y_sealevel)
    # Correct cell height if truncated
    cell_height = cell_height - truncated_cell_height
    return y_upper_left_cell, cell_height
end

""" Get x-coordinate of upper left corner of cell.

# Arguments
- `x_grid::Float64`: x-coordinate of grid node (meters)
- `cell_width::Float64`: Width of cell (meters)

# Returns
- `x_upper_left::Float64`: x-coordinate of upper left corner of cell (meters)
"""
function get_x_upper_left(x_grid::Float64, cell_width::Float64)::Float64
    return x_grid - cell_width / 2.0
end

""" Get y-coordinate of upper left corner of cell.

The cell is truncated by the sealevel.

# Arguments
- `y_grid::Float64`: y-coordinate of grid node (meters)
- `cell_height::Float64`: Height of cell (meters)
- `y_sealevel::Float64`: y-coordinate of sealevel (meters)

# Returns
- `y_upper_left::Float64`: y-coordinate of upper left corner of cell (meters)
- `truncated_cell_height::Float64`: Height of truncated portion (meters)
"""
function get_y_upper_left(
    y_grid::Float64, 
    cell_height::Float64, 
    y_sealevel::Float64
)::Tuple{Float64, Float64}
    y_upper_left = y_grid - cell_height / 2.0

    if y_upper_left < y_sealevel
        truncated_cell_height = y_sealevel - y_upper_left
    else
        truncated_cell_height = 0.0
    end

    y_upper_left = max(y_upper_left, y_sealevel)
    return y_upper_left, truncated_cell_height
end

""" Calculate horizontal cell geometry.

# Arguments
- `j::Int`: Index of grid node in x-direction
- `xstp::Vector{Float64}`: x-spacing of grid nodes (meters)
- `x_grid::Float64`: x-coordinate of grid node (meters)

# Returns
- `x_upper_left_cell::Float64`: x-coordinate of upper left corner of cell (meters)
- `cell_width::Float64`: Width of cell (meters)
"""
function calculate_horizontal_cell_geometry(
    j::Int,
    xstp::Vector{Float64},
    x_grid::Float64
)::Tuple{Float64, Float64}
    cell_width = get_cell_size_x(xstp, j)
    x_upper_left_cell = get_x_upper_left(x_grid, cell_width)
    return x_upper_left_cell, cell_width
end

""" Get cell size in x-direction.

# Arguments
- `xstp::Vector{Float64}`: x-spacing of grid nodes (meters)
- `j::Int`: Index of grid node in x-direction

# Returns
- `cell_size_x::Float64`: Cell size in x-direction (meters)
"""
function get_cell_size_x(xstp::Vector{Float64}, j::Int)::Float64
    nstep = length(xstp)
    if j == 1
        cell_size_x = xstp[j]
    elseif j < nstep
        cell_size_x = (xstp[j-1] + xstp[j]) / 2.0
    else
        cell_size_x = xstp[nstep]
    end
    return cell_size_x
end

""" Get cell size in y-direction.

# Arguments
- `ystp::Vector{Float64}`: y-spacing of grid nodes (meters)
- `i::Int`: Index of grid node in y-direction

# Returns
- `cell_size_y::Float64`: Cell size in y-direction (meters)
"""
function get_cell_size_y(ystp::Vector{Float64}, i::Int)::Float64
    nstep = length(ystp)
    if i == 1
        cell_size_y = ystp[i]
    elseif i < nstep
        cell_size_y = (ystp[i-1] + ystp[i]) / 2.0
    else
        cell_size_y = ystp[nstep]
    end
    return cell_size_y
end

end # module Gravity
