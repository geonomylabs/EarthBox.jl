module DepthDependentSideVelocity

import ..SolveForVelocity: solve_for_velocity
import ..VelocityFactors: calculate_velocity_factors

"""
    solve_for_veloc_bc_depth_dependent(
        y_linear::Float64,
        velocity_extension::Float64,
        ysize::Float64,
        xsize::Float64,
        sticky_thickness_left::Float64,
        sticky_thickness_right::Float64
    )::Tuple{Float64, Float64}

Calculate conservative top and bottom orthogonal velocity.

This function applies to cases where extension velocity applied to side
boundaries varies in a linear manner with depth below the lithosphere.
Extension velocity within the lithosphere is constant along side
boundaries.

# Arguments
- `y_linear::Float64`: Depth at which extension velocity becomes linear
- `velocity_extension::Float64`: Extension velocity
- `ysize::Float64`: Depth of model (m)
- `xsize::Float64`: Width of model (m)
- `sticky_thickness_left::Float64`: Thickness of sticky air on left side (m)
- `sticky_thickness_right::Float64`: Thickness of sticky air on right side (m)

# Returns
- Tuple containing:
  - `velocity_y_upper::Float64`:
    - Conservative vertical velocity (m/s) along upper boundary of extension model
  - `velocity_y_lower::Float64`:
    - Conservative vertical velocity (m/s) along lower boundary of extension model
"""
function solve_for_veloc_bc_depth_dependent(
    y_linear::Float64,
    velocity_extension::Float64,
    ysize::Float64,
    xsize::Float64,
    sticky_thickness_left::Float64,
    sticky_thickness_right::Float64
)::Tuple{Float64, Float64}
    vavg, f1_factor, f2_factor = depth_dependent_extension_velocity(
        y_linear, velocity_extension/2.0, ysize, xsize,
        sticky_thickness_left, sticky_thickness_right
    )

    velocity_extension_right = vavg/2.0
    velocity_extension_left = -vavg/2.0
    if abs(velocity_extension_left) > 0.0 && abs(velocity_extension_right) > 0.0
        velocity_y_upper, velocity_y_lower = solve_for_velocity(
            f1_factor, f2_factor, velocity_extension_right)
    else
        velocity_y_upper = 0.0
        velocity_y_lower = 0.0
    end
    return (velocity_y_upper, velocity_y_lower)
end

"""
    depth_dependent_extension_velocity(
        y_linear::Float64, 
        velocity_extension_half_rate::Float64, 
        ysize::Float64, 
        xsize::Float64, 
        sticky_thickness_left::Float64, 
        sticky_thickness_right::Float64
    )::Tuple{Float64, Float64, Float64}

Calculate parameters for depth-dependent extension model.

# Arguments
- `y_linear::Float64`: Depth at which extension velocity becomes linear
- `velocity_extension_half_rate::Float64`: Extension velocity at half rate (m/s)
- `ysize::Float64`: Depth of model (m)
- `xsize::Float64`: Width of model (m)
- `sticky_thickness_left::Float64`: Thickness of sticky air on left side (m)
- `sticky_thickness_right::Float64`: Thickness of sticky air on right side (m)

# Returns
- Tuple containing:
  - `vavg::Float64`: Average extension velocity (m/s)
  - `f1_factor::Float64`: Factor f1
  - `f2_factor::Float64`: Factor f2
"""
function depth_dependent_extension_velocity(
    y_linear::Float64,
    velocity_extension_half_rate::Float64,
    ysize::Float64,
    xsize::Float64,
    sticky_thickness_left::Float64,
    sticky_thickness_right::Float64
)::Tuple{Float64, Float64, Float64}
    velocity_average_half_rate = (
        (
            y_linear*velocity_extension_half_rate
            + (ysize - y_linear)*velocity_extension_half_rate*0.5
        )/ysize
    )
    f1_factor, f2_factor = calculate_velocity_factors(
        velocity_average_half_rate, xsize, ysize,
        sticky_thickness_left, sticky_thickness_right
    )
    return (velocity_average_half_rate, f1_factor, f2_factor)
end

end # module 