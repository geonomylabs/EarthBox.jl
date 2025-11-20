module ExtensionDepthDependentWithInflow

import EarthBox.ModelDataContainer: ModelData
import ..UpperAndLowerVelocity: calculate_conservative_upper_velocity,
                              calculate_conservative_lower_velocity
import ..PrintInflowOutflow: print_inflow_velocity, print_average_side_velocity
import ..GetInputs: get_extension_inputs
import ..SetParameters: set_velocity_upper_and_lower
import ..AverageVelocity: calculate_average_depth_dependent_velocity
import ..ConservationCheck: check_conservation_depth_dependent

"""
    calc_velocity!(model::ModelData; use_asymmetric_inflow_outflow::Bool=false)

Calculate velocity at top and bottom boundaries of extension models.

Boundary velocity is calculated to conserve mass based on extension velocity
along side boundaries.

# Arguments
- `y_linear`: Depth (m) at which extension velocity becomes linear.

# Returns
- Tuple containing `y_linear_left` and `y_linear_right`
"""
function calc_velocity!(
    model::ModelData;
    use_asymmetric_inflow_outflow::Bool=false
)::Tuple{Float64,Float64}
    debug = false # Set to True to print debug information

    extension_bc_inputs = get_extension_inputs(
        model, use_asymmetric_inflow_outflow=use_asymmetric_inflow_outflow)

    xsize = extension_bc_inputs.xsize
    ysize = extension_bc_inputs.ysize
    sticky_thickness_left = extension_bc_inputs.sticky_thickness_left
    sticky_thickness_right = extension_bc_inputs.sticky_thickness_right
    velocity_left = extension_bc_inputs.velocity_left
    velocity_right = extension_bc_inputs.velocity_right

    plate_thickness = model.bcs.parameters.velocity.plate_thickness.value

    y_linear_left = plate_thickness - sticky_thickness_left
    velocity_avg_left = calculate_average_depth_dependent_velocity(
        y_linear_left, velocity_left, ysize
    )

    y_linear_right = plate_thickness - sticky_thickness_right
    velocity_avg_right = calculate_average_depth_dependent_velocity(
        y_linear_right, velocity_right, ysize
    )

    if debug
        print_average_side_velocity(velocity_avg_left, velocity_avg_right)
    end

    velocity_y_upper = calculate_conservative_upper_velocity(
        xsize, velocity_left, velocity_right,
        sticky_thickness_left, sticky_thickness_right
    )

    velocity_y_lower = calculate_conservative_lower_velocity(
        ysize, xsize,
        velocity_avg_left, velocity_avg_right,
        velocity_y_upper
    )

    check_conservation_depth_dependent(
        xsize, ysize,
        y_linear_left,
        y_linear_right,
        velocity_avg_left, velocity_avg_right,
        velocity_y_upper, velocity_y_lower,
        velocity_left, velocity_right,
        sticky_thickness_left, sticky_thickness_right
    )

    if debug
        print_inflow_velocity(velocity_y_upper, velocity_y_lower)
    end

    set_velocity_upper_and_lower(model, velocity_y_upper, velocity_y_lower)

    return (y_linear_left, y_linear_right)
end

end # module 