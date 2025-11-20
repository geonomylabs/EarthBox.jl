module ContractionInflowAlongBottom

import EarthBox.ModelDataContainer: ModelData
import ..StickyThickness: get_sticky_thickness
import ..UpperAndLowerVelocity: calculate_conservative_upper_and_lower_velocity
import ..PrintInflowOutflow: print_inflow_velocity

"""
    calc_velocity!(model::ModelData; use_asymmetric_inflow_outflow::Bool=false)

Calculate velocity at top and bottom boundaries of contraction model.

Boundary velocity is calculated to conserve mass based contraction velocity along
side boundaries.

# Updated Parameters
- `model.bcs.parameters.velocity`
  - `vyu.value`: Conservative vertical velocity (m/s) along upper boundary of
    extension model.
  - `vyl.value`: Conservative vertical velocity (m/s) along lower boundary of
    extension model.
"""
function calc_velocity!(
    model::ModelData;
    use_asymmetric_inflow_outflow::Bool=false
)::Nothing
    debug = true

    full_velocity_contraction = abs(
        model.bcs.parameters.velocity.full_velocity_contraction.value
    )

    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value

    sticky_thickness_left, sticky_thickness_right = get_sticky_thickness(model)

    if !use_asymmetric_inflow_outflow
        velocity_left = full_velocity_contraction / 2.0
        velocity_right = -full_velocity_contraction / 2.0
    else
        velocity_left = 0.0
        velocity_right = -full_velocity_contraction
    end

    if full_velocity_contraction > 0.0
        (velocity_y_upper, velocity_y_lower) = 
            calculate_conservative_upper_and_lower_velocity(
                ysize,
                xsize,
                velocity_left,
                velocity_right,
                sticky_thickness_left,
                sticky_thickness_right
            )
    else
        velocity_y_upper = 0.0
        velocity_y_lower = 0.0
    end

    if debug
        print_inflow_velocity(velocity_y_upper, velocity_y_lower)
    end

    model.bcs.parameters.velocity.vyu.value = velocity_y_upper
    model.bcs.parameters.velocity.vyl.value = velocity_y_lower
    return nothing
end

end # module 