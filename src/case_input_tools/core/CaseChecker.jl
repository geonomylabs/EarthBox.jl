module CaseChecker

import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.PrintFuncs: print_info
import EarthBox.PrintFuncs: print_unit_conversion_info_and_values
import EarthBox: UnitConversion
import ..CaseTypes: CaseType
import ..CaseTypes: CaseCollectionType
import ..CaseTypes: CaseParameter

""" Check for required parameters and units.

# Arguments
- `required_inputs::Vector{Union{ParameterFloat, ParameterInt, ParameterStr}}`: List of required input parameters
- `case_inputs::CaseType`: Case inputs

# Throws
- `ArgumentError`: If required parameters are missing or have incorrect units
"""
function check_required_parameters(
    required_inputs::Vector{Union{ParameterFloat, ParameterInt, ParameterStr}},
    case_inputs::CaseType
)::Nothing
    for required_parameter in required_inputs
        if !haskey(case_inputs, required_parameter.name)
            required_names = [input.name for input in required_inputs]
            throw(ArgumentError(
                "Required parameters $(required_parameter.name) not found in case inputs. " *
                "Required parameters are: $(required_names)"
            ))
        end
        case_parameter_units = case_inputs[required_parameter.name].units
        if required_parameter.units != case_parameter_units
            throw(ArgumentError(
                "Parameter $(required_parameter.name) has units $(case_parameter_units) " *
                "but required units are $(required_parameter.units)"
            ))
        end
    end
end

""" 

    convert_parameters_to_standard_units!(
        parameters::Dict{String, CaseParameter},
        master_dict::Dict{String, Union{ParameterFloat, ParameterInt, ParameterStr}}
    )::Nothing

Update the units of parameters stored in the input dictionary called parameters to 
the standard units defined in the master_dict dictionary.

# Arguments
- `parameters::Dict{String, CaseParameter}`
    - Dictionary of parameters where the key is a parameter name and the value is 
       a CaseParameter struct containing the parameter value and units. Note this 
       function updates the units so they are consistent with standard units. 
- `master_dict::Dict{String, Union{ParameterFloat, ParameterInt, ParameterStr}}`
    - Dictionary of parameter information where the key is a parameter name and 
       the value is an Union{ParameterFloat, ParameterInt, ParameterStr} struct 
       containing the parameter type, default value, and units.
"""
function convert_parameters_to_standard_units!(
    parameters::Dict{String, CaseParameter},
    master_dict::Dict{String, Union{ParameterFloat, ParameterInt, ParameterStr}}
)::Nothing
    print_info("Converting case parameters to standard units")
    unit_conversion_data = UnitConversion.UnitConversionData()
    for (parameter_name, parameter_obj) in parameters
        unit_start_init = parameter_obj.units  # parameter units
        unit_end_init = master_dict[parameter_name].units  # standard units
        (
            conversion_func, unit_start, unit_end
        ) = UnitConversion.get_conversion_func(
            unit_start_init, unit_end_init, unit_conversion_data, parameter_name)
        value_start = parameter_obj.value
        value_end = conversion_func(value_start)
        parameters[parameter_name] = CaseParameter(value_end, unit_end)
        if unit_start != unit_end
            print_unit_conversion_info_and_values(
                parameter_name, value_start, value_end, unit_start, unit_end, level=2)
        end
    end
end

end # module