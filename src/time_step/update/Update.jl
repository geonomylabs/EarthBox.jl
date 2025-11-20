module Update

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData

""" Update main time step using marker displacement and velocity field.

This function calculates the new time step based on the maximum marker
displacement and velocity field and updates the model time step referred
to as timestep stored in the model data container.

# Updated Parameters
- `model.timestep.parameters.main_time_loop.timestep.value`: Main model timestep in seconds
"""
function update_model_time_step_using_displacement_limit(
    model::ModelData
)::Bool
    advection = model.markers.parameters.advection
    iuse_local_adaptive_time_stepping = advection.iuse_local_adaptive_time_stepping.value

    if iuse_local_adaptive_time_stepping == 0
        timestep_new, timestep_update = 
            calculate_new_time_step_using_average_grid_spacing(model)
    else
        print_info("Using local adaptive time stepping.", level=1)
        timestep_new, timestep_update = 
            calculate_timestep_using_local_staggered_grid_nodes(model)
    end

    update_timestep_in_domain(model, timestep_new)

    return timestep_update
end

function update_timestep_in_domain(
    model::ModelData,
    timestep_new::Float64
)::Nothing
    model.timestep.parameters.main_time_loop.timestep.value = timestep_new
    return nothing
end

function calculate_new_time_step_using_average_grid_spacing(
    model::ModelData
)::Tuple{Float64, Bool}
    vxmax, vymax = get_maximum_velocity_components(model)

    grid_geom = model.grids.parameters.geometry
    xstpavg = grid_geom.xstpavg.value
    ystpavg = grid_geom.ystpavg.value

    marker_cell_displ_max = model.markers.parameters.advection.marker_cell_displ_max.value
    timestep_old = model.timestep.parameters.main_time_loop.timestep.value

    timestep_new = calculate_new_timestep(
        timestep_old, marker_cell_displ_max, xstpavg, ystpavg, vxmax, vymax)

    timestep_update = get_timestep_update_boolean(timestep_old, timestep_new)

    return timestep_new, timestep_update
end

function get_timestep_update_boolean(
    timestep_old::Float64,
    timestep_new::Float64
)::Bool
    return timestep_new != timestep_old
end

function get_maximum_velocity_components(
    model::ModelData
)::Tuple{Float64, Float64}
    vel_arrays = model.stokes_continuity.arrays.staggered_grid_velocity
    vxmax = calculate_maximum_velocity_magnitude(vel_arrays.vx1.array)
    vymax = calculate_maximum_velocity_magnitude(vel_arrays.vy1.array)
    return vxmax, vymax
end

function calculate_maximum_velocity_magnitude(
    velocity_array::Array{Float64, 2}
)::Float64
    return max(abs(maximum(velocity_array)), abs(minimum(velocity_array)))
end

function calculate_new_timestep(
    timestep_old::Float64,
    marker_cell_displ_max::Float64,
    xstpavg::Float64,
    ystpavg::Float64,
    vxmax::Float64,
    vymax::Float64
)::Float64
    timestep_new = calculate_new_time_step_for_cell(
        timestep_old, marker_cell_displ_max, xstpavg, vxmax)
    timestep_new = calculate_new_time_step_for_cell(
        timestep_new, marker_cell_displ_max, ystpavg, vymax)
    return timestep_new
end

function calculate_timestep_using_local_staggered_grid_nodes(
    model::ModelData
)::Tuple{Float64, Bool}
    timestep_old = model.timestep.parameters.main_time_loop.timestep.value
    timestep_new = timestep_old

    marker_cell_displ_max = 
        model.markers.parameters.advection.marker_cell_displ_max.value

    vel_arrays = model.stokes_continuity.arrays.staggered_grid_velocity

    vx1 = vel_arrays.vx1.array
    xstp = model.grids.arrays.basic.xstp_b.array
    timestep_new = calculate_new_time_step_for_vx_staggered_grid(
        marker_cell_displ_max, vx1, xstp, timestep_new)

    vy1 = vel_arrays.vy1.array
    ystp = model.grids.arrays.basic.ystp_b.array
    timestep_new = calculate_new_time_step_for_vy_staggered_grid(
        marker_cell_displ_max, vy1, ystp, timestep_new)

    timestep_update = get_timestep_update_boolean(timestep_old, timestep_new)

    return timestep_new, timestep_update
end

function calculate_new_time_step_for_vx_staggered_grid(
    marker_cell_displ_max::Float64,
    vx1::Matrix{Float64},
    xstp::Vector{Float64},
    timestep::Float64
)::Float64
    for j in 1:size(vx1, 2)
        for i in 1:size(vx1, 1)
            cell_size = get_staggered_vx_cell_size(xstp, j)
            velocity_magnitude = abs(vx1[i, j])
            timestep = calculate_new_time_step_for_cell(
                timestep, marker_cell_displ_max,
                cell_size, velocity_magnitude
            )
        end
    end
    return timestep
end

@inline function get_staggered_vx_cell_size(
    xstp::Vector{Float64},
    j::Int
)::Float64
    jmax = length(xstp)
    if j == 1
        cell_size = xstp[1]
    elseif j == jmax + 1
        cell_size = xstp[end]
    else
        xstp_left = xstp[j-1]
        xstp_right = xstp[j]
        cell_size = min(xstp_left, xstp_right)
    end
    return cell_size
end

function calculate_new_time_step_for_vy_staggered_grid(
    marker_cell_displ_max::Float64,
    vy1::Matrix{Float64},
    ystp::Vector{Float64},
    timestep::Float64
)::Float64
    for j in 1:size(vy1, 2)
        for i in 1:size(vy1, 1)
            cell_size = get_staggered_vy_cell_size(ystp, i)
            velocity_magnitude = abs(vy1[i, j])
            timestep = calculate_new_time_step_for_cell(
                timestep, marker_cell_displ_max,
                cell_size, velocity_magnitude
            )
        end
    end
    return timestep
end

@inline function get_staggered_vy_cell_size(
    ystp::Vector{Float64},
    i::Int
)::Float64
    imax = length(ystp)
    if i == 1
        cell_size = ystp[1]
    elseif i == imax + 1
        cell_size = ystp[end]
    else
        ystep_up = ystp[i-1]
        ystep_down = ystp[i]
        cell_size = min(ystep_up, ystep_down)
    end
    return cell_size
end

@inline function calculate_new_time_step_for_cell(
    timestep::Float64,
    marker_cell_displ_max::Float64,
    cell_size::Float64,
    velocity_magnitude::Float64
)::Float64
    if velocity_magnitude > 0
        timestep_limit = calculate_time_step_limit(
            marker_cell_displ_max, cell_size, velocity_magnitude)
        timestep = min(timestep, timestep_limit)
    end
    return timestep
end

@inline function calculate_time_step_limit(
    marker_cell_displ_max::Float64,
    cell_size::Float64,
    velocity_magnitude::Float64
)::Float64
    return marker_cell_displ_max * cell_size / velocity_magnitude
end

end # module 