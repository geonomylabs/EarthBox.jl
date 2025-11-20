"""
Functions used to calculate the thermal structure of the oceanic lithosphere
"""
module HalfSpaceCooling

import SpecialFunctions: erf

"""
    calculate_temperature(
        depth_below_surface_meters::Float64,
        surface_temperature_kelvin::Float64,
        asthenosphere_temperature_kelvin::Float64,
        distance_from_ridge_meters::Float64,
        half_spreading_rate_m_s::Float64,
        thermal_diffusivity_m2_s::Float64
    )

Calculate temperature at a given depth below the seafloor.

This function uses the half-space cooling model to calculate the
temperature.

# Arguments
- `depth_below_surface_meters::Float64`: Depth below the seafloor in meters.
- `surface_temperature_kelvin::Float64`: Surface temperature in Kelvin.
- `asthenosphere_temperature_kelvin::Float64`: Asthenosphere temperature in Kelvin.
- `distance_from_ridge_meters::Float64`: Distance from the ridge in meters.
- `half_spreading_rate_m_s::Float64`: Half spreading rate in meters per second.
- `thermal_diffusivity_m2_s::Float64`: Thermal diffusivity in square meters per second.

# Returns
- `temperature_kelvin::Float64`: 
    - Temperature at the given depth below the seafloor in Kelvin.
"""
function calculate_temperature(
    depth_below_surface_meters::Float64,
    surface_temperature_kelvin::Float64,
    asthenosphere_temperature_kelvin::Float64,
    distance_from_ridge_meters::Float64,
    half_spreading_rate_m_s::Float64,
    thermal_diffusivity_m2_s::Float64
)::Float64
    temperature_kelvin = (
        surface_temperature_kelvin
        + (asthenosphere_temperature_kelvin - surface_temperature_kelvin)
        * erf(
            depth_below_surface_meters / 2.0
            / sqrt(
                thermal_diffusivity_m2_s
                * distance_from_ridge_meters / half_spreading_rate_m_s
            )
        )
    )
    return temperature_kelvin
end

end # module 