module FormatValueStrings

import EarthBox: ConversionFuncs
using Printf

"""
    get_value_string_for_files(param_name::String, value::Union{Float64,Int64}, units_str::String)::Tuple{String,String}

Define string for file output.
"""
function get_value_string_for_files(
    param_name::String,
    value::Union{Float64,Int64},
    units_str::String
)::Tuple{String,String}
    value, units_str = convert_to_input_file_units(param_name, value, units_str)
    if value isa Int64
        value_str = @sprintf("%12d", value)
    else
        if value >= 1e6 || value <= 1e-3
            value_str = @sprintf("%12.4e", value)
        else
            value_str = @sprintf("%12.4f", value)
        end
    end
    return value_str, units_str
end

"""
    convert_to_input_file_units(param_name::String, value::Union{Float64,Int64}, units_str::String)::Tuple{Union{Float64,Int64},String}

Convert from model units to input file units.
"""
function convert_to_input_file_units(
    param_name::String,
    value::Union{Float64,Int64},
    units_str::String
)::Tuple{Union{Float64,Int64},String}
    if units_str == "m/s"
        if param_name in ["erosion_rate", "sedimentation_rate"]
            value = ConversionFuncs.meters_per_seconds_to_mm_per_yr(value)
            units_str = "mm/yr"
        else
            value = ConversionFuncs.meters_per_seconds_to_cm_per_yr(value)
            units_str = "cm/yr"
        end
    elseif units_str == "K" && param_name != "max_temp_change"
        value = ConversionFuncs.kelvin_to_celsius(value)
        units_str = "C"
    elseif units_str == "s"
        value = ConversionFuncs.seconds_to_years(value)
        units_str = "yr"
    end
    return value, units_str
end

"""
    get_value_string_for_terminal(value::Union{Float64,Int64,String})::String

Define string for terminal output.
"""
function get_value_string_for_terminal(
    value::Union{Float64,Int64,String}
)::String
    if value isa Int64
        value_str = @sprintf("%-12d", value)
    elseif value isa Float64
        if 1e-4 < value < 1e6
            value_str = @sprintf("%-12.4f", value)
        else
            value_str = @sprintf("%-12.4e", value)
        end
    else
        value_str = string(value)
    end
    return value_str
end

end # module 