module Sediment

"""
    melting_model_parameters(pressure_pascals::Float64)::Tuple{Float64,Float64,Float64}

Calculate liquidus and solidus temperatures and define latent heat.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins.
- `temperature_solidus::Float64`: Solidus temperature in Kelvins.
- `latent_heat::Float64`: Latent heat in J/kg.
"""
function melting_model_parameters(pressure_pascals::Float64)::Tuple{Float64,Float64,Float64}
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, MPa
    pressure_mpa = pressure_pascals * 1e-6
    if pressure_mpa < 1200.0
        temperature_solidus = (
            889.0 +
            17900.0 / (pressure_mpa + 54.0) +
            20200.0 / (pressure_mpa + 54.0)^2.0
        )
    else
        temperature_solidus = 831.0 + 0.06 * pressure_mpa
    end
    temperature_liquidus = 1262.0 + 0.09 * pressure_mpa
    latent_heat = 300000.0
    return (temperature_liquidus, temperature_solidus, latent_heat)
end

"""
    liquidus_gerya2010(pressure_pascals::Float64)::Float64

Calculate liquidus temperature for sediment.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins.
"""
function liquidus_gerya2010(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, MPa
    pressure_mpa = pressure_pascals * 1e-6
    temperature_liquidus = 1262.0 + 0.09 * pressure_mpa
    return temperature_liquidus
end

"""
    solidus_gerya2010(pressure_pascals::Float64)::Float64

Calculate solidus temperature for sediment.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_solidus::Float64`: Solidus temperature in Kelvins.
"""
function solidus_gerya2010(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, MPa
    pressure_mpa = pressure_pascals * 1e-6
    if pressure_mpa < 1200.0
        temperature_solidus = (
            889.0 +
            17900.0 / (pressure_mpa + 54.0) +
            20200.0 / (pressure_mpa + 54.0)^2.0
        )
    else
        temperature_solidus = 831.0 + 0.06 * pressure_mpa
    end
    return temperature_solidus
end

end # module Sediment 