module ExtensionInflowAndOutflowAlongSides

import EarthBox.ModelDataContainer: ModelData
import ..UpperAndLowerVelocity: calculate_conservative_upper_velocity
import ..GetInputs: get_left_and_right_extension_velocities
import ..AverageVelocity: calculate_average_outflow_velocity
import ..SetParameters: set_velocity_upper, set_inflow_velocity
import ..OutflowGeometry: get_outflow_geometry, get_outflow_size
import ..InflowGeometry: get_inflow_geometry
import ..StickyThickness: get_sticky_thickness
import ..ConservationCheck: check_conservation_inflow_and_outflow_along_sides

"""
    calc_velocity!(model::ModelData; use_asymmetric_inflow_outflow::Bool=false)

Calculate velocity at top and side inflow boundaries of extension model.

Boundary velocity is calculated to conserve volume in the model domain.

# Updated Parameters
- `model.bcs.parameters.velocity`
  - `vyu.value`: Conservative vertical velocity (m/s) along upper boundary of
    extension model.
  - `velocity_inflow_left.value`: Inflow velocity (m/s) along left boundary of
    extension model.
  - `velocity_inflow_right.value`: Inflow velocity (m/s) along right boundary of
    extension model.
"""
function calc_velocity!(
    model::ModelData;
    use_asymmetric_inflow_outflow::Bool=false
)::Nothing
    (sticky_thickness_left, sticky_thickness_right) = get_sticky_thickness(model)

    velocity_left, velocity_right = get_left_and_right_extension_velocities(
        model, use_asymmetric_inflow_outflow=use_asymmetric_inflow_outflow
    )

    xsize = model.grids.parameters.geometry.xsize.value
    velocity_y_upper = calculate_conservative_upper_velocity(
        xsize, velocity_left, velocity_right,
        sticky_thickness_left, sticky_thickness_right
    )
    set_velocity_upper(model, velocity_y_upper)

    (velocity_avg_outflow_left, velocity_avg_outflow_right) = 
        calculate_average_outflow_velocities(model, velocity_left, velocity_right)

    outflow = calculate_outflow(
        model, velocity_avg_outflow_left, velocity_avg_outflow_right
    )

    (velocity_inflow_left, velocity_inflow_right) = calculate_inflow_velocities(
        model, outflow, velocity_y_upper
    )
    set_inflow_velocity(model, velocity_inflow_left, velocity_inflow_right)

    check_conservation_inflow_and_outflow_along_sides(
        model, velocity_y_upper,
        velocity_avg_outflow_left, velocity_avg_outflow_right,
        velocity_left, velocity_right
    )
end

"""
    calculate_average_outflow_velocities(model::ModelData, velocity_left::Float64,
                                       velocity_right::Float64)::Tuple{Float64,Float64}

Calculate average outflow velocities.
"""
function calculate_average_outflow_velocities(
    model::ModelData,
    velocity_left::Float64,
    velocity_right::Float64
)::Tuple{Float64,Float64}
    plate_thickness, smoothing_thickness = get_outflow_geometry(model)
    sticky_thickness_left, sticky_thickness_right = get_sticky_thickness(model)

    velocity_avg_outflow_left = calculate_average_outflow_velocity(
        sticky_thickness_left, plate_thickness, smoothing_thickness,
        velocity_left
    )

    velocity_avg_outflow_right = calculate_average_outflow_velocity(
        sticky_thickness_right, plate_thickness, smoothing_thickness,
        velocity_right
    )

    return (velocity_avg_outflow_left, velocity_avg_outflow_right)
end

"""
    calculate_outflow(model::ModelData, velocity_avg_outflow_left::Float64,
                     velocity_avg_outflow_right::Float64)::Float64

Calculate outflow along side boundaries.
"""
function calculate_outflow(
    model::ModelData,
    velocity_avg_outflow_left::Float64,
    velocity_avg_outflow_right::Float64
)::Float64
    outflow_size_left, outflow_size_right = get_outflow_size(model)
    outflow = (
        velocity_avg_outflow_right * outflow_size_right -
        velocity_avg_outflow_left * outflow_size_left
    )

    return outflow
end

"""
    calculate_inflow_velocities(model::ModelData, outflow::Float64,
                              velocity_y_upper::Float64)::Tuple{Float64,Float64}

Calculate inflow velocities.

# Arguments
- `outflow`: Outflow flux along side boundaries (left + right).
- `velocity_y_upper`: Conservative vertical velocity along upper boundary.
"""
function calculate_inflow_velocities(
    model::ModelData,
    outflow::Float64,
    velocity_y_upper::Float64
)::Tuple{Float64,Float64}
    _, smoothing_thickness = get_outflow_geometry(model)
    ysize_inflow_left, ysize_inflow_right = get_inflow_geometry(model)

    xsize = model.grids.parameters.geometry.xsize.value
    half_side_inflow = (outflow - velocity_y_upper * xsize) / 2.0

    ysize_inflow_left_constant = ysize_inflow_left - smoothing_thickness
    ysize_inflow_right_constant = ysize_inflow_right - smoothing_thickness

    velocity_inflow_left = (
        half_side_inflow /
        (smoothing_thickness * 0.5 + ysize_inflow_left_constant)
    )

    velocity_inflow_right = (
        -half_side_inflow /
        (smoothing_thickness * 0.5 + ysize_inflow_right_constant)
    )

    return (velocity_inflow_left, velocity_inflow_right)
end

end # module 