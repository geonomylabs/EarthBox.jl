module CaseBuilder

import EarthBox.ParameterRegistry: get_eb_parameters
import ..CaseTypes: CaseCollectionType, CaseParameter

""" 
    define_case_group!(
        case_inputs::CaseCollectionType;
        case_id_ini::Int,
        parameter_name::String,
        values::Vector{Float64},
        units::String,
        fixed_parameter_name::Union{String, Nothing}=nothing,
        fixed_value::Union{Float64, Nothing}=nothing,
    )::Int

Define case inputs for a given target parameter name and a list of values. 
An optional fixed key can be set to a fixed value for all cases.

# Arguments
- `case_inputs::CaseCollectionType`
    - Dictionary containing the case inputs for each case.
- `case_id_ini::Int`
    - Initial case ID to start building cases from.
- `parameter_name::String`
    - Parameter name in the case inputs dictionary to modify for each case
       using values from the `values` vector.
- `values::Vector{Float64}`
    - List of values to assign to the parameter name for each case.
- `units::String`
    - Units to assign to the parameter name for each case.
- `fixed_parameter_name::Union{String, Nothing}`
    - Parameter name key to set a fixed value for all cases.
- `fixed_value::Union{Float64, Nothing}`
    - Fixed value to assign to the fixed_parameter_name for all cases.
- `fixed_units::Union{String, Nothing}`
    - Units to assign to the fixed_parameter_name for all cases.

# Returns
- `case_id::Int`: The case ID of the last case built

"""
function define_case_group!(
    case_inputs::CaseCollectionType;
    case_id_ini::Int,
    parameter_name::String,
    values::Vector{Float64},
    units::String,
    fixed_parameter_name::Union{String, Nothing}=nothing,
    fixed_value::Union{Float64, Nothing}=nothing,
    fixed_units::Union{String, Nothing}=nothing
)::Int
    keys = get_eb_parameters()
    valid_keys = [getfield(keys, f).name for f in fieldnames(typeof(get_eb_parameters()))]
    
    if !(parameter_name in valid_keys)
        throw(ArgumentError(
            "Invalid target key: $parameter_name. Expected one of $valid_keys"
        ))
    end
    
    if fixed_parameter_name !== nothing
        if !(fixed_parameter_name in valid_keys)
            throw(ArgumentError(
                "Invalid fixed key: $fixed_parameter_name. Expected one of $valid_keys"
            ))
        end
        if fixed_units == nothing
            throw(ArgumentError("Fixed units cannot be nothing if fixed key is set."))
        end
    end
    
    if isempty(values)
        throw(ArgumentError("Values list cannot be empty."))
    end
    
    icount = 0
    for value in values
        case_id = case_id_ini + icount
        case = "case$case_id"
        case_inputs[case][parameter_name] = CaseParameter(value, units)
        if fixed_parameter_name !== nothing && fixed_value !== nothing
            case_inputs[case][fixed_parameter_name] = CaseParameter(fixed_value, fixed_units)
        end
        icount += 1
    end
    
    return case_id_ini + icount - 1
end

end # module 