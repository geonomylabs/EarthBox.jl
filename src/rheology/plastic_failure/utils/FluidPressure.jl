module FluidPressure

"""
    calculate_fluid_pressure(y_marker::Float64, y_sealevel::Float64, density_fluid::Float64=1000.0)::Float64

Calculate fluid pressure.

# Arguments
- `y_marker::Float64`: Y-coordinate of marker in m
- `y_sealevel::Float64`: Y-coordinate of sealevel in m
- `density_fluid::Float64`: Density of fluid in kg/m^3 (default: 1000.0)

# Returns
- `pressure_fluid::Float64`: Fluid pressure in Pa
"""
@inline function calculate_fluid_pressure(
    y_marker::Float64,
    y_sealevel::Float64,
    density_fluid::Float64=1000.0
)::Float64
    gravity = 9.81
    depth_below_sealevel = max(0.0, y_marker - y_sealevel)
    pressure_fluid = density_fluid * gravity * depth_below_sealevel
    return pressure_fluid
end

end # module 