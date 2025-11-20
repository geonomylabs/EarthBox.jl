module Gabbro

"""
    solidus_gerya2010_gabbro_glacier(pressure_pascals::Float64)::Float64

Modification of the dry gabbro solidus for the gabbro glacier model.

See Maclennan et al. (2004) figure 5, which shows a shift in the solidus
relative to the Gerya 2010 model of around 100C based on the MELTS model.

This shift in the solidus is associated with crystal fractionation in
shallow magma chambers.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_solidus::Float64`: Solidus temperature in Kelvins.
"""
function solidus_gerya2010_gabbro_glacier(pressure_pascals::Float64)::Float64
    temperature_solidus = solidus_gerya2010(pressure_pascals)
    temperature_solidus = temperature_solidus + 100.0
    return temperature_solidus
end

"""
    liquidus_gerya2010_gabbro_glacier(pressure_pascals::Float64)::Float64

Modification of the dry gabbro liquidus for the gabbro glacier model.

See Maclennan et al. (2004) figure 5, which shows a shift in the liquidus
relative to the Gerya 2010 model of around 300C based on the MELTS model.

This shift in the liquidus is associated with crystal fractionation in
shallow magma chambers.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins.
"""
function liquidus_gerya2010_gabbro_glacier(pressure_pascals::Float64)::Float64
    temperature_liquidus = liquidus_gerya2010(pressure_pascals)
    temperature_liquidus = temperature_liquidus + 300.0
    return temperature_liquidus
end

"""
    solidus_gerya2010(pressure_pascals::Float64)::Float64

Calculate the dry solidus temperature.

Dry gabbro solidus from Hess (1989) as summarized in Gerya 2010/2013.

This function applies to gabbroic magma and solidified gabbro.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_solidus::Float64`: Solidus temperature in Kelvins.
"""
function solidus_gerya2010(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, MPa
    pressure_mpa = pressure_pascals * 1e-6
    temperature_solidus = 1327.0 + 0.091 * pressure_mpa
    return temperature_solidus
end

"""
    liquidus_gerya2010(pressure_pascals::Float64)::Float64

Calculate liquidus temperature.

This function applies to gabbroic magma and solidified gabbro.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins.
"""
function liquidus_gerya2010(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, MPa
    pressure_mpa = pressure_pascals * 1e-6
    # Liquidus temperature
    # From Gerya 2013 (Dry) with reference to Hess (1989).
    temperature_liquidus = 1423.0 + 0.105 * pressure_mpa
    return temperature_liquidus
end

"""
    solidus_schubert2013(pressure_pascals::Float64)::Float64

Calculate the wet solidus temperature.

Wet Gabbro melting model from Schmidt and Poli (1998).

This function applies to gabbroic magma and solidified gabbro.

# References
Schmidt and Poli (1998) in "Experimentally based water budgets for dehydrating
slabs and consequences for arc magma generation", Earth and Planetary Science
Letters, 163, 361-379.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_solidus::Float64`: Solidus temperature in Kelvins.
"""
function solidus_schubert2013(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, MPa
    pressure_mpa = pressure_pascals * 1e-6

    # Schubert et al. (2013), Table 1 with reference to Hess (1989),
    # Schmidt and Poli (1998), Hirschmann (2000), Johannes (1985) and Poli
    # and Schmidt (2002).
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

end # module Gabbro 