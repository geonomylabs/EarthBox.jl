module Coordinates

import Random
import EarthBox.ModelDataContainer: ModelData
import ..Sticky: check_recycle_sticky
import ..Sticky: get_sticky_material_ids
import ...BoundaryVelocity.InflowRamp: calculate_ramp_factor

"""
    RecycledCoordinates

Structure holding new recycled coordinates in meters for sticky and subsurface markers.
"""
struct RecycledCoordinates
    newx_air::Vector{Float64}
    newy_air::Vector{Float64}
    newx_subsurface::Vector{Float64}
    newy_subsurface::Vector{Float64}
end

"""
    MarkerRecycleArrays

Structure including recycled_indices, an array storing marker indices that will be recycled
and recycle_flags, an array of flags indicating if markers are recycled in air (0) or
recycled subsurface material (1).
"""
struct MarkerRecycleArrays
    recycle_indices::Vector{Int64}
    recycle_flags::Vector{Int64}
end

"""
    InflowOutflowInputs

Inputs used to calculate recycled marker coordinates for cases where vertical boundaries
have inflow and/or outflow regions.
"""
Base.@kwdef struct InflowOutflowInputs
    nrecycle::Int64 = 0
    nrecycle_start::Int64 = 0
    y_height::Float64 = 0.0
    y_initial::Float64 = 0.0
    x_initial::Float64 = 0.0
    y_smooth_top::Float64 = 0.0
    y_smooth_bottom::Float64 = 0.0
    displacement_x::Float64 = 0.0
end

""" Update recycled marker coordinates.

# Arguments
- `model`: Model data container
- `recycled_coordinates`: New coordinates for recycled markers
- `outside_indices`: Indices of markers outside the domain

# Returns
- `MarkerRecycleArrays`: Structure containing recycled indices and flags
"""
function reset_recycled_marker_coordinates(
    model::ModelData,
    recycled_coordinates::RecycledCoordinates,
    outside_indices::Vector{Int64}
)::MarkerRecycleArrays
    newx_air = recycled_coordinates.newx_air
    newy_air = recycled_coordinates.newy_air
    newx_subsurface = recycled_coordinates.newx_subsurface
    newy_subsurface = recycled_coordinates.newy_subsurface

    nrecycle_air = length(newx_air)
    nrecycle_subsurface = length(newx_subsurface)

    marker_arrays = model.markers.arrays
    location = marker_arrays.location
    marker_x = location.marker_x.array
    marker_y = location.marker_y.array
    marker_matid = marker_arrays.material.marker_matid.array

    matids_sticky = get_sticky_material_ids(model)

    nrecycle_total = nrecycle_air + nrecycle_subsurface
    recycle_indices = zeros(Int64, nrecycle_total)
    recycle_flags = zeros(Int64, nrecycle_total)

    noutside = length(outside_indices)

    icount_air = 1
    icount_subsurface = 1
    for ioutside in 1:noutside
        imarker = outside_indices[ioutside]
        recycle_sticky = check_recycle_sticky(
            marker_matid[imarker], matids_sticky, nrecycle_air)
        if recycle_sticky
            icount_air = update_marker_coordinates!(
                icount_air, newx_air, newy_air,
                imarker, marker_x, marker_y,
            )
            recycle_flag = 0 # Recycled from sticky air/water
        elseif nrecycle_subsurface > 0
            icount_subsurface = update_marker_coordinates!(
                icount_subsurface, newx_subsurface, newy_subsurface,
                imarker, marker_x, marker_y,
            )
            recycle_flag = 1 # Recycled from subsurface
        end
        recycle_indices[ioutside] = imarker
        recycle_flags[ioutside] = recycle_flag
    end

    return MarkerRecycleArrays(
        recycle_indices, 
        recycle_flags
    )
end

function update_marker_coordinates!(
    icount::Int64,
    newx::Vector{Float64},
    newy::Vector{Float64},
    imarker::Int64,
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
)::Int64
    marker_x[imarker] = newx[icount]
    marker_y[imarker] = newy[icount]
    return icount + 1
end

function calc_recycled_x_coordinates_horizontal_boundary(
    nrecycle::Int64,
    x_initial::Float64,
    xwidth::Float64
)::Vector{Float64}
    newx = Vector{Float64}(undef, nrecycle)
    for j in 1:nrecycle
        random_num = rand()
        newx[j] = x_initial + random_num * xwidth
    end
    return newx
end

function calc_recycled_y_coordinates_horizontal_boundary(
    nrecycle::Int64,
    y_initial::Float64,
    displacement_y::Float64
)::Vector{Float64}
    newy = Vector{Float64}(undef, nrecycle)
    for j in 1:nrecycle
        random_num = rand()
        newy[j] = y_initial + random_num * displacement_y
    end
    return newy
end

function update_recycled_coordinates_vertical_inflow_outflow_boundary!(
    inflow_outflow_inputs::InflowOutflowInputs,
    newy::Vector{Float64},
    newx::Vector{Float64}
)::Nothing
    nrecycle = inflow_outflow_inputs.nrecycle
    nrecycle_start = inflow_outflow_inputs.nrecycle_start
    y_height = inflow_outflow_inputs.y_height
    y_initial = inflow_outflow_inputs.y_initial
    x_initial = inflow_outflow_inputs.x_initial
    y_smooth_top = inflow_outflow_inputs.y_smooth_top
    y_smooth_bottom = inflow_outflow_inputs.y_smooth_bottom
    displacement_x = inflow_outflow_inputs.displacement_x
    calculate_recycled_y_coordinates_vertical_boundary!(
        nrecycle, nrecycle_start, y_height, y_initial, newy)
    calculate_recycled_x_coordinates_vertical_inflow_outflow_boundary!(
        nrecycle, nrecycle_start, x_initial, y_smooth_top, y_smooth_bottom, 
        displacement_x, newy, newx)
    return nothing
end

"""
    calculate_recycled_x_coordinates_vertical_inflow_outflow_boundary!(
        nrecycle::Int64,
        nrecycle_start::Int64,
        x_initial::Float64,
        y_smooth_top::Float64,
        y_smooth_bottom::Float64,
        displacement_x::Float64,
        newy::Vector{Float64},
        newx::Vector{Float64}
    )::Nothing

Calculate recycled x-locations for markers on vertical boundary with inflow and outflow.

# Updated arrays
- 'newx::Vector{Float64}' (nrecycle_right)
    - New recycled x-locations for subsurface markers.
"""
function calculate_recycled_x_coordinates_vertical_inflow_outflow_boundary!(
    nrecycle::Int64,
    nrecycle_start::Int64,
    x_initial::Float64,
    y_smooth_top::Float64,
    y_smooth_bottom::Float64,
    displacement_x::Float64,
    newy::Vector{Float64},
    newx::Vector{Float64}
)::Nothing
    Threads.@threads for j in 1:nrecycle
        random_num = rand()
        mm = j + nrecycle_start
        ramp_factor = calculate_ramp_factor(newy[mm], y_smooth_top, y_smooth_bottom)
        newx[mm] = x_initial + random_num * displacement_x * ramp_factor
    end
    return nothing
end

function update_recycled_coordinates_vertical_boundary!(
    inflow_outflow_inputs::InflowOutflowInputs,
    newy::Vector{Float64},
    newx::Vector{Float64}
)::Nothing
    nrecycle = inflow_outflow_inputs.nrecycle
    nrecycle_start = inflow_outflow_inputs.nrecycle_start
    y_height = inflow_outflow_inputs.y_height
    y_initial = inflow_outflow_inputs.y_initial
    x_initial = inflow_outflow_inputs.x_initial
    displacement_x = inflow_outflow_inputs.displacement_x
    calculate_recycled_y_coordinates_vertical_boundary!(
        nrecycle, nrecycle_start, y_height, y_initial, newy)

    calculate_recycled_x_coordinates_vertical_boundary!(
        nrecycle, nrecycle_start, x_initial, displacement_x, newx)

    return nothing
end

"""
    calculate_recycled_y_coordinates_vertical_boundary!(
        nrecycle::Int64,
        nrecycle_start::Int64,
        y_height::Float64,
        y_initial::Float64,
        newy::Vector{Float64}
    )::Nothing

Calculate recycled y-locations for subsurface markers on vertical boundary.

# Updated arrays
- 'newy::Vector{Float64}' (nrecycle)
    - New recycled y-locations for subsurface markers.
"""
function calculate_recycled_y_coordinates_vertical_boundary!(
    nrecycle::Int64,
    nrecycle_start::Int64,
    y_height::Float64,
    y_initial::Float64,
    newy::Vector{Float64}
)::Nothing
    Threads.@threads for j in 1:nrecycle
        random_num = rand()
        mm = j + nrecycle_start
        newy[mm] = y_initial + random_num * y_height
    end
    return nothing
end

"""
    calculate_recycled_x_coordinates_vertical_boundary!(
        nrecycle::Int64,
        nrecycle_start::Int64,
        x_initial::Float64,
        displacement_x::Float64,
        newx::Vector{Float64}
    )::Nothing

Calculate recycled air x-coordinates for vertical boundaries.

# Updated arrays
- 'newx::Vector{Float64}' (nrecycle_left)
    - New recycled x-locations markers.
"""
function calculate_recycled_x_coordinates_vertical_boundary!(
    nrecycle::Int64,
    nrecycle_start::Int64,
    x_initial::Float64,
    displacement_x::Float64,
    newx::Vector{Float64}
)::Nothing
    Threads.@threads for j in 1:nrecycle
        random_num = rand()
        mm = j + nrecycle_start
        newx[mm] = x_initial + random_num * displacement_x
    end
    return nothing
end

function get_random_numbers_for_recycled_markers(
    nrecycle::Int64
)::Vector{Float64}
    random_numbers = Vector{Float64}(undef, nrecycle)
    Threads.@threads for j in 1:nrecycle
        random_numbers[j] = rand()
    end
    return random_numbers
end

end # module Coordinates 