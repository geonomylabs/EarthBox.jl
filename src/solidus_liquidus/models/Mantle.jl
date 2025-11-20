module Mantle

"""
    solidus_gerya2010(pressure_pascals::Float64)::Float64

Calculate solidus temperature for mantle.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_solidus::Float64`: Solidus temperature in Kelvins.
"""
function solidus_gerya2010(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, MPa
    pressure_mpa = pressure_pascals * 1e-6

    # (Gerya 2010, Table 17.2)
    # Lithospheric mantle (dry)
    if pressure_mpa < 10000.0
        temperature_solidus = (
            1393.811 +
            0.132899 * pressure_mpa -
            0.000005104 * pressure_mpa^2.0
        )
    else
        temperature_solidus = (
            2212.4 +
            0.030819 * (pressure_mpa - 10000.0)
        )
    end
    return temperature_solidus
end

"""
    liquidus_gerya2010(pressure_pascals::Float64)::Float64

Calculate liquidus temperature for mantle.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins.
"""
function liquidus_gerya2010(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, MPa
    pressure_mpa = pressure_pascals * 1e-6
    temperature_liquidus = 2073.0 + 0.114 * pressure_mpa
    return temperature_liquidus
end

"""
    solidus_katz2003(pressure_pascals::Float64)::Float64

Calculate solidus temperature using Katz 2003 model.

# References
Katz, R. F., M. Spiegelman, and C. H. Langmuir, A new parameterization of
hydrous mantle melting, Geochem. Geophys. Geosyst.,4(9), 1073,
doi:10.1029/2002GC000433, 2003.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_solidus::Float64`: Solidus temperature in Kelvins.
"""
function solidus_katz2003(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, GPa
    pressure_gpa = pressure_pascals * 1e-9

    # Dry Peridotite from Katz (2003)
    if pressure_gpa <= 8.0
        # The Katz (2003) model is only applicable to ~ 8 GPa
        temperature_solidus_celsius = solidus_katz2003_anhydrous_less_than_8gpa(pressure_gpa)
    else
        # use a linearization similar to Gerya (201)
        temperature_solidus_celsius_at_8gpa = solidus_katz2003_anhydrous_less_than_8gpa(8.0)
        temperature_solidus_celsius_at_20gpa = 2250.0
        temperature_solidus_celsius = linear_above_8pga(
            pressure_gpa,
            temperature_solidus_celsius_at_8gpa,
            temperature_solidus_celsius_at_20gpa
        )
    end
    temperature_solidus = temperature_solidus_celsius + 273.0
    return temperature_solidus
end

"""
    solidus_katz2003_anhydrous_less_than_8gpa(pressure_gpa::Float64)::Float64

Calculate solidus temperature for pressure less than 8 GPa using Katz 2003 model.

# Arguments
- `pressure_gpa::Float64`: Pressure in GPa.

# Returns
- `temperature_solidus_celsius::Float64`: Solidus temperature in Celsius.
"""
function solidus_katz2003_anhydrous_less_than_8gpa(pressure_gpa::Float64)::Float64
    temperature_solidus_celsius = (
        -5.1 * pressure_gpa^2.0 +
        132.9 * pressure_gpa +
        1085.7
    )
    return temperature_solidus_celsius
end

"""
    liquidus_katz2003(pressure_pascals::Float64)::Float64

Calculate liquidus temperature using Katz 2003 model.

# References
Katz, R. F., M. Spiegelman, and C. H. Langmuir, A new parameterization of
hydrous mantle melting, Geochem. Geophys. Geosyst.,4(9), 1073,
doi:10.1029/2002GC000433, 2003.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals.

# Returns
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins.
"""
function liquidus_katz2003(pressure_pascals::Float64)::Float64
    pressure_pascals = max(pressure_pascals, 0.0)
    # Pressure, GPa
    pressure_gpa = pressure_pascals * 1e-9

    if pressure_gpa <= 8.0
        temperature_liquidus_celsius = liquidus_katz2003_anhydrous_less_than_8gpa(pressure_gpa)
    else
        temperature_liquidus_celsius_at_8gpa = liquidus_katz2003_anhydrous_less_than_8gpa(8.0)
        temperature_liquidus_celsius_at_20gpa = 2250.0 + 100.0
        temperature_liquidus_celsius = linear_above_8pga(
            pressure_gpa,
            temperature_liquidus_celsius_at_8gpa,
            temperature_liquidus_celsius_at_20gpa
        )
    end
    temperature_liquidus = temperature_liquidus_celsius + 273.0
    return temperature_liquidus
end

"""
    liquidus_katz2003_anhydrous_less_than_8gpa(pressure_gpa::Float64)::Float64

Calculate liquidus temperature for pressure less than 8 GPa using Katz 2003 model.

# Arguments
- `pressure_gpa::Float64`: Pressure in GPa.

# Returns
- `temperature_liquidus_celsius::Float64`: Liquidus temperature in Celsius.
"""
function liquidus_katz2003_anhydrous_less_than_8gpa(pressure_gpa::Float64)::Float64
    temperature_liquidus_celsius = (
        -2.0 * pressure_gpa^2.0 +
        45.0 * pressure_gpa +
        1780.0
    )
    return temperature_liquidus_celsius
end

"""
    linear_above_8pga(pressure_gpa::Float64, temperature_celsius_at_8gpa::Float64, temperature_celsius_at_20gpa::Float64)::Float64

Linearize solidus/liquidus temperature above 8 GPa.

# Arguments
- `pressure_gpa::Float64`: Pressure in GPa.
- `temperature_celsius_at_8gpa::Float64`: Temperature at 8 GPa in Celsius.
- `temperature_celsius_at_20gpa::Float64`: Temperature at 20 GPa in Celsius.

# Returns
- `temperature_celsius::Float64`: Linearized temperature in Celsius.
"""
function linear_above_8pga(
    pressure_gpa::Float64,
    temperature_celsius_at_8gpa::Float64,
    temperature_celsius_at_20gpa::Float64
)::Float64
    temperature_celsius = (
        temperature_celsius_at_8gpa +
        (temperature_celsius_at_20gpa - temperature_celsius_at_8gpa) *
        (pressure_gpa - 8.0) / 12.0
    )
    return temperature_celsius
end

end # module Mantle 