module PressureReset

"""
    calculate_reset_pressure(
        y_marker::Float64;
        density::Float64=3300.0,
        gravity::Float64=9.81
    ) -> Float64

Calculate pressure reset value.

# Arguments
- `y_marker::Float64`: y-location of marker (m)
- `density::Float64=3300.0`: Density of material (kg/m^3)
- `gravity::Float64=9.81`: Acceleration due to gravity (m/s^2)

# Returns
- `Float64`: Calculated pressure reset value
"""
function calculate_reset_pressure(
    y_marker::Float64;
    density::Float64=3300.0,
    gravity::Float64=9.81
)::Float64
    return density * gravity * y_marker
end

end # module