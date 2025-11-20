module UniformSideVelocity

import ..SolveForVelocity: solve_for_velocity
import ..VelocityFactors: calculate_velocity_factors

"""
    solve_for_veloc_bc_uniform_side_velocity(
        velocity_magnitude_half_rate::Float64, ysize::Float64, xsize::Float64, 
        sticky_thickness_left::Float64, sticky_thickness_right::Float64; 
        contraction::Bool=false
    )::Tuple{Float64, Float64}

Calculate conservative top and bottom orthogonal velocity.

This function applies to cases with uniform extension or contraction
velocity along side boundaries.

# Arguments
- `velocity_magnitude_half_rate::Float64`: 
    - Half rate extension or contraction velocity (m/s) along side boundaries
- `ysize::Float64`: 
    - Vertical dimension of model (m)
- `xsize::Float64`: 
    - Horizontal dimension of model (m)
- `sticky_thickness_left::Float64`: 
    - Thickness of sticky air layer along left boundary (m)
- `sticky_thickness_right::Float64`: 
    - Thickness of sticky air layer along right boundary (m)
- `contraction::Bool=false`: 
    - Flag to indicate contraction model. If false extension will be assumed

# Returns
- Tuple containing:
  - `velocity_y_upper::Float64`:
    - Conservative vertical velocity (m/s) along upper boundary of extension model
  - `velocity_y_lower::Float64`:
    - Conservative vertical velocity (m/s) along lower boundary of extension model
"""
function solve_for_veloc_bc_uniform_side_velocity(
    velocity_magnitude_half_rate::Float64,
    ysize::Float64,
    xsize::Float64,
    sticky_thickness_left::Float64,
    sticky_thickness_right::Float64;
    contraction::Bool=false
)::Tuple{Float64, Float64}
    velocity_magnitude_half_rate = abs(velocity_magnitude_half_rate)
    f1_factor, f2_factor = calculate_velocity_factors(
        velocity_magnitude_half_rate, xsize, ysize,
        sticky_thickness_left, sticky_thickness_right
    )
    velocity_y_upper, velocity_y_lower = solve_for_velocity(
        f1_factor, f2_factor, velocity_magnitude_half_rate)
    if contraction
        velocity_y_upper = -velocity_y_upper
        velocity_y_lower = -velocity_y_lower
    end
    return (velocity_y_upper, velocity_y_lower)
end

end # module UniformSideVelocity 