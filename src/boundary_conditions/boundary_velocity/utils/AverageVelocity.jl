module AverageVelocity

"""
    calculate_average_depth_dependent_velocity(
        y_linear::Float64, plate_velocity::Float64, ysize::Float64)::Float64

Calculate average depth-dependent velocity along side boundary.

# Arguments
- `y_linear::Float64`: Depth at which extension velocity becomes linear (m)
- `plate_velocity::Float64`: Uniform velocity in plate above linear depth-dependent zone (m/s)
- `ysize::Float64`: Depth of model (m)

# Returns
- `velocity_boundary_average::Float64`: Average velocity along boundary (m/s)
"""
function calculate_average_depth_dependent_velocity(
    y_linear::Float64,
    plate_velocity::Float64,
    ysize::Float64
)::Float64
    velocity_boundary_average = (
        (
            y_linear*plate_velocity
            + (ysize - y_linear)*plate_velocity*0.5
        )/ysize
    )
    return velocity_boundary_average
end

"""
    calculate_average_sub_plate_depth_dependent_velocity(plate_velocity::Float64)::Float64

Calculate average depth-dependent velocity along side boundary.

# Arguments
- `plate_velocity::Float64`: Uniform velocity in plate above linear depth-dependent zone (m/s)

# Returns
- `velocity_below_plate_average::Float64`: Average velocity below plate (m/s)
"""
function calculate_average_sub_plate_depth_dependent_velocity(
    plate_velocity::Float64
)::Float64
    velocity_below_plate_average = plate_velocity*0.5
    return velocity_below_plate_average
end

"""
    calculate_average_outflow_velocity(sticky_thickness::Float64, plate_thickness::Float64, smoothing_thickness::Float64, plate_velocity::Float64)::Float64

Calculate average outflow velocity.

# Arguments
- `sticky_thickness::Float64`: Thickness of sticky zone (m)
- `plate_thickness::Float64`: Thickness of plate (m)
- `smoothing_thickness::Float64`: Thickness of smoothing zone (m)
- `plate_velocity::Float64`: Uniform velocity in plate above linear depth-dependent zone (m/s)

# Returns
- `velocity_outflow_average::Float64`: Average velocity along outflow boundary (m/s)
"""
function calculate_average_outflow_velocity(
    sticky_thickness::Float64,
    plate_thickness::Float64,
    smoothing_thickness::Float64,
    plate_velocity::Float64
)::Float64
    part1 = sticky_thickness*plate_velocity
    part2 = plate_thickness*plate_velocity
    part3 = smoothing_thickness*plate_velocity*0.5

    outflow_thickness = (
        sticky_thickness + plate_thickness + smoothing_thickness)

    velocity_outflow_average = (part1 + part2 + part3)/outflow_thickness

    return velocity_outflow_average
end

end # module AverageVelocity 