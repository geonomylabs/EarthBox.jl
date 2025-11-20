module ContinentalCrustLower

""" Calculate solidus temperature for lower continental crust.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_solidus::Float64`: Solidus temperature in Kelvins.
"""
function solidus_gerya2010(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, MPa
    pressure_mpa = pressure_pascals * 1e-6
    if pressure_mpa < 1600.0
        temperature_solidus = (
            973.0 -
            70400.0 / (pressure_mpa + 354.0) +
            77800000.0 / (pressure_mpa + 354.0)^2.0
        )
    else
        temperature_solidus = (
            935.0 +
            0.0035 * pressure_mpa +
            0.0000062 * pressure_mpa^2.0
        )
    end
    return temperature_solidus
end

""" Calculate liquidus temperature for lower continental crust.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins.
"""
function liquidus_gerya2010(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, MPa
    pressure_mpa = pressure_pascals * 1e-6
    temperature_liquidus = 1423.0 + 0.105 * pressure_mpa
    return temperature_liquidus
end

end # module ContinentalCrustLower 