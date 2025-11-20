module VelocityFactors

"""
    calculate_velocity_factors(
        velocity_magnitude_half_rate::Float64,
        xsize::Float64,
        ysize::Float64,
        sticky_thickness_left::Float64, 
        sticky_thickness_right::Float64
    )::Tuple{Float64, Float64}

Calculate factors for solving for inflow velocity.

# Arguments
- `velocity_magnitude_half_rate`: Half of the velocity magnitude rate
- `xsize`: Size in x direction
- `ysize`: Size in y direction
- `sticky_thickness_left`: Sticky thickness on the left side
- `sticky_thickness_right`: Sticky thickness on the right side

# Returns
- Tuple containing f1_factor and f2_factor
"""
function calculate_velocity_factors(
    velocity_magnitude_half_rate::Float64,
    xsize::Float64,
    ysize::Float64,
    sticky_thickness_left::Float64,
    sticky_thickness_right::Float64
)::Tuple{Float64, Float64}
    f1_factor = 2.0 * velocity_magnitude_half_rate * ysize / xsize
    constant1 = (sticky_thickness_left + sticky_thickness_right) / 2.0 / ysize
    f2_factor = (constant1 - 1.0) / constant1
    return (f1_factor, f2_factor)
end

end # module VelocityFactors 