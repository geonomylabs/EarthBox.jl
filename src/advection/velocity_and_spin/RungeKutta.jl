module RungeKutta

import StaticArrays: MVector

""" Calculate x coordinate for the next Runge-Kutta cycle.

No need to calc new x_marker_current for rk = 3 since this is the last iteration.
"""
@inline function calculate_x_for_next_runge_kutta_cycle(
    rk_order::Int,
    timestep::Float64,
    x_marker_initial::Float64,
    velocity_x::Float64
)::Float64
    if rk_order < 4
        if rk_order < 3  # rk_order 1 and 2
            return x_marker_initial + timestep / 2.0 * velocity_x
        else  # rk_order 3
            return x_marker_initial + timestep * velocity_x
        end
    end
    return x_marker_initial  # Default return for rk_order >= 3
end

""" Calculate y coordinate for the next Runge-Kutta cycle.

No need to calc new y_marker_current for rk = 3 since this is the last iteration.
"""
@inline function calculate_y_for_next_runge_kutta_cycle(
    rk_order::Int,
    timestep::Float64,
    y_marker_initial::Float64,
    velocity_y::Float64
)::Float64
    if rk_order < 4
        if rk_order < 3  # rk_order 1 and 2
            return y_marker_initial + timestep / 2.0 * velocity_y
        else  # rk_order 3
            return y_marker_initial + timestep * velocity_y
        end
    end
    return y_marker_initial  # Default return for rk_order >= 3
end

""" Calculate final velocity or spin using 4-th order Runge-Kutta.

# Arguments
- `rkarray::Vector{Float64}`: Array containing velocity components (m/s) or spin (1/s)
    interpolated at 4 points in the velocity field.

# Returns
- `Float64`: Interpolated velocity component or spin using 4-th order Runge-Kutta.
"""
@inline function apply_4th_order_runge_kutta(rkarray::Vector{Float64})::Float64
    @inbounds return (rkarray[1] + 2.0 * rkarray[2] + 2.0 * rkarray[3] + rkarray[4]) / 6.0
end

@inline function apply_4th_order_runge_kutta(rkarray::MVector{4, Float64})::Float64
    @inbounds return (rkarray[1] + 2.0 * rkarray[2] + 2.0 * rkarray[3] + rkarray[4]) / 6.0
end

@inline function apply_4th_order_runge_kutta(
    vx_r1::Float64, 
    vx_r2::Float64, 
    vx_r3::Float64, 
    vx_r4::Float64
)::Float64
    @inbounds return (vx_r1 + 2.0 * vx_r2 + 2.0 * vx_r3 + vx_r4) / 6.0
end

end # module 