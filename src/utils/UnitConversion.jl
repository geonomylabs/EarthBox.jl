module UnitConversion

import EarthBox.ConversionFuncs

const UnitFuncsType = Dict{String, Dict{String, Function}}

const UnitKeysType = Dict{String, Vector{String}}

"""
    UnitConversionData

Container for unit conversion functions and registry.

# Fields
- `unit_funcs::UnitFuncsType`: Dictionary mapping from source unit to dictionary of target units and conversion functions
- `unit_registry::UnitKeysType`: Dictionary mapping unit categories to valid unit strings

# Constructor
    UnitConversionData()::UnitConversionData

# Returns
- `UnitConversionData`: New UnitConversionData with initialized conversion functions and registry
"""
struct UnitConversionData
    unit_funcs::UnitFuncsType
    unit_registry::UnitKeysType
end

function UnitConversionData()
    unit_funcs = get_unit_funcs()
    unit_registry = get_unit_registry()
    return UnitConversionData(unit_funcs, unit_registry)
end

"""
    get_unit_funcs()::UnitFuncsType

Get the unit conversion functions dictionary where the key is the standard unit 
name and the value is a dictionary of conversion functions for converting to 
other units.

# Standard Units and Available Conversions

| Standard Unit | Available Conversions |
|---------------|-----------------------|
| `cm/yr`       | `m/s`, `cm/hr`        |
| `mm/yr`       | `m/s`                 |
| `m/yr`        | `m/s`                 |
| `cm/hr`       | `m/s`                 |
| `C`           | `K`                   |
| `K`           | `C`                   |
| `yr`          | `s`, `Myr`            |
| `s`           | `yr`, `Myr`           |
| `Myr`         | `yr`, `s`             |
| `m`           | `km`, `cm`            |
| `km`          | `m`, `cm`             |
| `Pa.s`        | `log10(Pa.s)`         |
| `1/s`         | `log10(1/s)`          |
| `Pa`          | `MPa`, `GPa`          |
| `m^2/yr`      | `m^2/s`               |
| `percent`     | `fraction`            |
| `deltaK`      | `deltaC`              |
| `deltaC`      | `deltaK`              |

# Returns
- `UnitFuncsType`: Dictionary of unit conversion functions.
"""
function get_unit_funcs()::UnitFuncsType
    unit_funcs = Dict{String, Dict{String, Function}}(
        "cm/yr" => Dict(
            "m/s" => ConversionFuncs.cm_per_yr_to_meters_per_seconds,
            "cm/hr" => cm_per_yr_to_cm_per_hour
        ),
        "mm/yr" => Dict(
            "m/s" => ConversionFuncs.mm_per_yr_to_meters_per_seconds
        ),
        "m/yr" => Dict(
            "m/s" => ConversionFuncs.meters_per_year_to_meters_per_seconds
        ),
        "cm/hr" => Dict(
            "m/s" => cm_per_hr_to_meters_per_second
        ),
        "C" => Dict(
            "K" => ConversionFuncs.celsius_to_kelvin
        ),
        "K" => Dict(
            "C" => ConversionFuncs.kelvin_to_celsius
        ),
        "yr" => Dict(
            "s" => ConversionFuncs.years_to_seconds
        ),
        "s" => Dict(
            "yr" => ConversionFuncs.seconds_to_years,
            "Myr" => seconds_to_millions_of_years
        ),
        "Myr" => Dict(
            "yr" => millions_of_years_to_years,
            "s" => ConversionFuncs.millions_of_years_to_seconds
        ),
        "m" => Dict(
            "km" => meters_to_kilometers,
            "cm" => meters_to_centimeters
        ),
        "km" => Dict(
            "m" => kilometers_to_meters,
            "cm" => kilometers_to_centimeters
        ),
        "Pa.s" => Dict(
            "log10(Pa.s)" => log10
        ),
        "1/s" => Dict(
            "log10(1/s)" => log10
        ),
        "Pa" => Dict(
            "MPa" => x -> x/1e6,
            "GPa" => x -> x/1e9
        ),
        "m^2/yr" => Dict(
            "m^2/s" => ConversionFuncs.meters_squared_per_year_to_meters_squared_per_second
        ),
        "percent" => Dict(
            "fraction" => x -> x/100.0
        ),
        "deltaK" => Dict(
            "deltaC" => x -> x
        ),
        "deltaC" => Dict(
            "deltaK" => x -> x
        )
    )
    return unit_funcs
end

"""
    get_unit_registry()::UnitKeysType

Get dictionary of standard unit names and valid alternative names.
Any valid alternative will be converted to the standard unit when a
parameter is loaded into the model.

# Standard Units and Valid Alternative Names

| Standard Unit   | Valid Alternative Names                         |
|-----------------|-------------------------------------------------|
| `m`             | `meters`                                        |
| `cm`            | `centimeters`                                   |
| `km`            | `kilometers`                                    |
| `m^3`           | `cubic meters`                                  |
| `s`             | `seconds`                                       |
| `yr`            | `years`                                         |
| `Myr`           | `millions of years`                             |
| `Ma`            | `megaannum`                                     |
| `kg`            | `kilograms`                                     |
| `kg/m^3`        | `kilograms per cubic meters`                    |
| `m/s`           | `meters/second`                                 |
| `m/yr`          | `meters per year`                               |
| `cm/yr`         | `centimeters per year`                          |
| `cm/hr`         | `centimeters per hour`                          |
| `mm/yr`         | `millimeters per year`                          |
| `1/s`           | `per second`                                    |
| `log10(1/s)`    | `log10(strain rate)`                            |
| `K`             | `kelvin`                                        |
| `C`             | `Celsius`                                       |
| `1/K`           | `per kelvin`                                    |
| `m/s/s`         | `m/s^2`                                         |
| `m^2/s`         | `meters squared per second`                     |
| `m^2/yr`        | `meters squared per year`                       |
| `Pa`            | `pascals`                                       |
| `MPa`           | `megapascals`                                   |
| `GPa`           | `gigapascals`                                   |
| `1/Pa`          | `per Pascal`                                    |
| `Pa.s`          | `pascal seconds`                                |
| `log10(Pa.s)`   | `log10(pascal seconds)`                         |
| `J`             | `joules`                                        |
| `W`             | `watts`                                         |
| `W/m`           | `watts per meter`                               |
| `K/m`           | `Kelvins per meter`                             |
| `K/km`          | `Kelvins per kilometer`                         |
| `C/m`           | `Celsius per meter`                             |
| `J/kg/K`        | `joules per kilogram per kelvin`                |
| `W/m/K`         | `watts per meter per kelvin`                    |
| `W/m^3`         | `watts per cubic meter`                         |
| `percent`       | `%`                                             |
| `fraction`      | `fraction`, `frac`                              |
| `degrees`       | `Degrees`, `Degree`                             |
| `1/s/MPa^n`     | `per second per megapascal to the nth power`    |
| `s^-m1*MPa^-m2` | `s to the -m1 power times MPa to the -m2 power` |
| `kJ/mol`        | `kilojoules per mole`                           |
| `J/MPa/mol`     | `joules per megapascal per mole`                |
| `J/kg`          | `joules per kilogram`                           |
| `deltaK`        | `delta Kelvin`                                  |
| `deltaC`        | `delta Celsius`                                 |

# Returns
- `UnitKeysType`: Dictionary of standard unit names and valid alternative names.

"""
function get_unit_registry()::UnitKeysType
    unit_registry = Dict{String, Vector{String}}(
        "m" => ["meters"],
        "cm" => ["centimeters"],
        "km" => ["kilometers"],
        "m^3" => ["cubic meters"],
        "s" => ["seconds"],
        "yr" => ["years"],
        "Myr" => ["millions of years"],
        "Ma" => ["megaannum"],
        "kg" => ["kilograms"],
        "kg/m^3" => ["kilograms per cubic meters"],
        "m/s" => ["meters/second"],
        "m/yr" => ["meters per year"],
        "cm/yr" => ["centimeters per year"],
        "cm/hr" => ["centimeters per hour"],
        "mm/yr" => ["millimeters per year"],
        "1/s" => ["per second"],
        "log10(1/s)" => ["log10(strain rate)"],
        "K" => ["kelvin"],
        "C" => ["Celsius"],
        "1/K" => ["per kelvin"],
        "m/s/s" => ["m/s^2"],
        "m^2/s" => ["meters squared per second"],
        "m^2/yr" => ["meters squared per year"],
        "Pa" => ["pascals"],
        "MPa" => ["megapascals"],
        "GPa" => ["gigapascals"],
        "1/Pa" => ["per Pascal"],
        "Pa.s" => ["pascal seconds"],
        "log10(Pa.s)" => ["log10(pascal seconds)", "log10(Pas)"],
        "J" => ["joules"],
        "W" => ["watts"],
        "W/m" => ["watts per meter"],
        "K/m" => ["Kelvins per meter"],
        "K/km" => ["Kelvins per kilometer"],
        "C/m" => ["Celsius per meter"],
        "J/kg/K" => ["joules per kilogram per kelvin"],
        "W/m/K" => ["watts per meter per kelvin"],
        "W/m^3" => ["watts per cubic meter"],
        "percent" => ["%"],
        "fraction" => ["fraction", "frac"],
        "degrees" => ["Degrees", "Degree"],
        "1/s/MPa^n" => [],
        "s^-m1*MPa^-m2" => [],
        "kJ/mol" => ["kilojoules per mole"],
        "J/MPa/mol" => ["joules per megapascal per mole"],
        "J/kg" => ["joules per kilogram"],
        "deltaK" => ["delta Kelvin"],
        "deltaC" => ["delta Celsius"],
        "None" => ["none"]
    )
    return unit_registry
end

function get_standard_unit(unit::String, unit_registry::UnitKeysType)::String
    standard_names = collect(keys(unit_registry))
    if unit in standard_names
        return unit
    end
    
    for (standard_unit, alternative_names) in unit_registry
        if unit in alternative_names
            return standard_unit
        end
    end
    
    throw(ArgumentError("The unit name $unit is not registered."))
end

function get_conversion_func(
    unit_start::String,
    unit_end::String,
    unit_conversion_data::UnitConversionData,
    parameter_name::String
)::Tuple{Function, String, String}
    unit_start = get_standard_unit(unit_start, unit_conversion_data.unit_registry)
    unit_end = get_standard_unit(unit_end, unit_conversion_data.unit_registry)
    
    if unit_start == unit_end
        return (no_conversion, unit_start, unit_end)
    end
    
    if !haskey(unit_conversion_data.unit_funcs, unit_start)
        throw(ArgumentError(
            "The starting unit $unit_start for parameter $parameter_name is not present in the unit_funcs dictionary associated with $unit_end."
        ))
    end
    
    if !haskey(unit_conversion_data.unit_funcs[unit_start], unit_end)
        throw(ArgumentError(
            "The end unit $unit_end for parameter $parameter_name is not present in the unit_funcs dictionary associated with $unit_start."
        ))
    end
    
    return (unit_conversion_data.unit_funcs[unit_start][unit_end], unit_start, unit_end)
end

function cm_per_yr_to_cm_per_hour(velocity_cm_yr::Float64)::Float64
    return velocity_cm_yr/(365.0*24.0)
end

function cm_per_hr_to_meters_per_second(velocity_cm_hr::Float64)::Float64
    velocity_cm_yr = velocity_cm_hr*365.0*24.0
    return ConversionFuncs.cm_per_yr_to_meters_per_seconds(velocity_cm_yr)
end

function seconds_to_millions_of_years(seconds::Float64)::Float64
    years = ConversionFuncs.seconds_to_years(seconds)
    return years/1e6
end

function millions_of_years_to_years(millions_of_years::Float64)::Float64
    return millions_of_years*1e6
end

function meters_to_kilometers(meters::Float64)::Float64
    return meters/1000.0
end

function meters_to_centimeters(meters::Float64)::Float64
    return meters*100.0
end

function kilometers_to_meters(kilometers::Float64)::Float64
    return kilometers*1000.0
end

function kilometers_to_centimeters(kilometers::Float64)::Float64
    return kilometers*1000.0*100.0
end

function log10(value::Float64)::Float64
    return value > 0.0 ? Base.log10(value) : 0.0
end

function no_conversion(value::Union{Float64, Int64, String, Symbol})::Union{Float64, Int64, String, Symbol}
    return value
end

end # module 