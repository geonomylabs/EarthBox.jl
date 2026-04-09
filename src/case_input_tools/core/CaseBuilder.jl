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
        fixed_units::Union{String, Nothing}=nothing,
        fixed_parameter_name2::Union{String, Nothing}=nothing,
        fixed_value2::Union{Float64, Nothing}=nothing,
        fixed_units2::Union{String, Nothing}=nothing,
        fixed_parameter_name3::Union{String, Nothing}=nothing,
        fixed_value3::Union{Float64, Nothing}=nothing,
        fixed_units3::Union{String, Nothing}=nothing,
        fixed_parameter_name4::Union{String, Nothing}=nothing,
        fixed_value4::Union{Float64, Nothing}=nothing,
        fixed_units4::Union{String, Nothing}=nothing,
    )::Int

Define case inputs for a given target parameter name and a list of values. 
Optional fixed keys (up to four) can be set to fixed values for all cases.

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
- `fixed_parameter_name`, `fixed_parameter_name2`, …, `fixed_parameter_name4` (`Union{String, Nothing}`)
    - Parameter name keys to set to fixed values for all cases.
- `fixed_value`, `fixed_value2`, …, `fixed_value4` (`Union{Float64, Nothing}`)
    - Fixed values for the corresponding fixed parameter names.
- `fixed_units`, `fixed_units2`, …, `fixed_units4` (`Union{String, Nothing}`)
    - Units for the corresponding fixed parameter names.

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
    fixed_units::Union{String, Nothing}=nothing,
    fixed_parameter_name2::Union{String, Nothing}=nothing,
    fixed_value2::Union{Float64, Nothing}=nothing,
    fixed_units2::Union{String, Nothing}=nothing,
    fixed_parameter_name3::Union{String, Nothing}=nothing,
    fixed_value3::Union{Float64, Nothing}=nothing,
    fixed_units3::Union{String, Nothing}=nothing,
    fixed_parameter_name4::Union{String, Nothing}=nothing,
    fixed_value4::Union{Float64, Nothing}=nothing,
    fixed_units4::Union{String, Nothing}=nothing,
)::Int
    keys = get_eb_parameters()
    valid_keys = [getfield(keys, f).name for f in fieldnames(typeof(get_eb_parameters()))]
    
    if !(parameter_name in valid_keys)
        throw(ArgumentError(
            "Invalid target key: $parameter_name. Expected one of $valid_keys"
        ))
    end

    fixed_triples = (
        (fixed_parameter_name, fixed_value, fixed_units),
        (fixed_parameter_name2, fixed_value2, fixed_units2),
        (fixed_parameter_name3, fixed_value3, fixed_units3),
        (fixed_parameter_name4, fixed_value4, fixed_units4),
    )
    for (fpn, _fv, fu) in fixed_triples
        if fpn !== nothing
            if !(fpn in valid_keys)
                throw(ArgumentError(
                    "Invalid fixed key: $fpn. Expected one of $valid_keys"
                ))
            end
            if fu === nothing
                throw(ArgumentError("Fixed units cannot be nothing if fixed key is set."))
            end
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
        for (fpn, fv, fu) in fixed_triples
            if fpn !== nothing && fv !== nothing
                case_inputs[case][fpn] = CaseParameter(fv, fu)
            end
        end
        icount += 1
    end
    
    return case_id_ini + icount - 1
end

end # module 