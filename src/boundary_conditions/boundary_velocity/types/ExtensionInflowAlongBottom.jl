module ExtensionInflowAlongBottom

import EarthBox.ModelDataContainer: ModelData
import ..UpperAndLowerVelocity: calculate_conservative_upper_and_lower_velocity
import ..PrintInflowOutflow: print_inflow_velocity
import ..GetInputs: get_extension_inputs
import ..SetParameters: set_velocity_upper_and_lower
import ..ConservationCheck: check_conservation_constant

"""
    calc_velocity!(model::ModelData; use_asymmetric_inflow_outflow::Bool=false)

Calculate velocity at top and bottom boundaries of extension model.

Boundary velocity is calculated to conserve mass based on extension velocity
along side boundaries.
"""
function calc_velocity!(
    model::ModelData;
    use_asymmetric_inflow_outflow::Bool=false
)::Nothing
    debug = false

    extension_bc_inputs = get_extension_inputs(
        model, use_asymmetric_inflow_outflow=use_asymmetric_inflow_outflow)

    xsize = extension_bc_inputs.xsize
    ysize = extension_bc_inputs.ysize
    sticky_thickness_left = extension_bc_inputs.sticky_thickness_left
    sticky_thickness_right = extension_bc_inputs.sticky_thickness_right
    velocity_left = extension_bc_inputs.velocity_left
    velocity_right = extension_bc_inputs.velocity_right

    (
        velocity_y_upper, 
        velocity_y_lower
    ) = calculate_conservative_upper_and_lower_velocity(
            ysize, xsize, velocity_left, velocity_right,
            sticky_thickness_left, sticky_thickness_right
        )

    check_conservation_constant(
        xsize, ysize, velocity_y_upper, velocity_y_lower,
        velocity_left, velocity_right,
        sticky_thickness_left, sticky_thickness_right
    )

    if debug
        print_inflow_velocity(velocity_y_upper, velocity_y_lower)
    end

    set_velocity_upper_and_lower(model, velocity_y_upper, velocity_y_lower)
    return nothing
end

end # module 