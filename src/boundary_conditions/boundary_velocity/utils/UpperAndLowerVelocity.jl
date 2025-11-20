module UpperAndLowerVelocity

"""
    calculate_conservative_upper_and_lower_velocity(
        ysize::Float64, xsize::Float64, velocity_left::Float64, 
        velocity_right::Float64, sticky_thickness_left::Float64, 
        sticky_thickness_right::Float64
    )::Tuple{Float64, Float64}

Calculate conservative top and bottom orthogonal velocity.

This function applies to cases with uniform extension or contraction
velocity along side boundaries.

# Arguments
- `ysize::Float64`: Vertical dimension of model (m)
- `xsize::Float64`: Horizontal dimension of model (m)
- `velocity_left::Float64`: Velocity along left boundary (m/s)
- `velocity_right::Float64`: Velocity along right boundary (m/s)
- `sticky_thickness_left::Float64`: Thickness of sticky air layer along left boundary (m)
- `sticky_thickness_right::Float64`: Thickness of sticky air layer along right boundary (m)

# Returns
- Tuple containing:
  - `velocity_y_upper::Float64`:
    - Conservative vertical velocity (m/s) along upper boundary of extension model
  - `velocity_y_lower::Float64`:
    - Conservative vertical velocity (m/s) along lower boundary of extension model
"""
function calculate_conservative_upper_and_lower_velocity(
    ysize::Float64,
    xsize::Float64,
    velocity_left::Float64,
    velocity_right::Float64,
    sticky_thickness_left::Float64,
    sticky_thickness_right::Float64
)::Tuple{Float64, Float64}
    velocity_y_upper = calculate_conservative_upper_velocity(
        xsize, velocity_left, velocity_right,
        sticky_thickness_left, sticky_thickness_right
    )
    velocity_y_lower = calculate_conservative_lower_velocity(
        ysize, xsize, velocity_left, velocity_right, velocity_y_upper
    )
    return (velocity_y_upper, velocity_y_lower)
end

"""
    calculate_conservative_upper_velocity(
        xsize::Float64, velocity_left::Float64, velocity_right::Float64, 
        sticky_thickness_left::Float64, sticky_thickness_right::Float64
    )::Float64

Calculate conservative top and bottom orthogonal velocity.

This function applies to cases with uniform extension or contraction
velocity along side boundaries.

# Arguments
- `xsize::Float64`: Horizontal dimension of model (m)
- `velocity_left::Float64`: Velocity along the sticky left boundary (m/s)
- `velocity_right::Float64`: Velocity along sticky right boundary (m/s)
- `sticky_thickness_left::Float64`: Thickness of sticky air layer along left boundary (m)
- `sticky_thickness_right::Float64`: Thickness of sticky air layer along right boundary (m)

# Returns
- `velocity_y_upper::Float64`:
    - Conservative vertical velocity (m/s) along upper boundary of extension model
"""
function calculate_conservative_upper_velocity(
    xsize::Float64,
    velocity_left::Float64,
    velocity_right::Float64,
    sticky_thickness_left::Float64,
    sticky_thickness_right::Float64
)::Float64
    full_velocity_extension = abs(velocity_left) + abs(velocity_right)
    if full_velocity_extension > 0.0
        velocity_y_upper = (
            (
                velocity_right*sticky_thickness_right
                - velocity_left*sticky_thickness_left
            )/xsize
        )
    else
        velocity_y_upper = 0.0
    end
    return velocity_y_upper
end

"""
    calculate_conservative_lower_velocity(
        ysize::Float64, xsize::Float64, velocity_left::Float64, 
        velocity_right::Float64, velocity_y_upper::Float64
    )::Float64

Calculate conservative top and bottom orthogonal velocity.

This function applies to cases with uniform extension or contraction
velocity along side boundaries.

# Arguments
- `ysize::Float64`: Vertical dimension of model (m)
- `xsize::Float64`: Horizontal dimension of model (m)
- `velocity_left::Float64`: Average velocity along left boundary (m/s)
- `velocity_right::Float64`: Average velocity along right boundary (m/s)
- `velocity_y_upper::Float64`: Conservative vertical velocity (m/s) along upper boundary of extension model

# Returns
- `velocity_y_lower::Float64`:
    - Conservative vertical velocity (m/s) along lower boundary of extension model
"""
function calculate_conservative_lower_velocity(
    ysize::Float64,
    xsize::Float64,
    velocity_left::Float64,
    velocity_right::Float64,
    velocity_y_upper::Float64
)::Float64
    full_velocity_extension = abs(velocity_left) + abs(velocity_right)
    if full_velocity_extension > 0.0
        velocity_y_lower = (
            velocity_y_upper
            - (velocity_right - velocity_left)*ysize/xsize
        )
    else
        velocity_y_lower = 0.0
    end
    return velocity_y_lower
end

end # module UpperAndLowerVelocity 