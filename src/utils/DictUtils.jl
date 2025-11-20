module DictUtils

import EarthBox.EarthBoxDtypes: InputDictType

"""
    update_initialization_dict(input_dict::Dict{String, Union{Float64, Int}},
        required_params::Dict{String, Function}, use_lists::Bool = true)::Dict{String, Any}

Updates an initialization dictionary by converting values to their required types. 
Checks that all required parameters are present in the input dictionary and applies 
type conversion functions to each value.

# Arguments
- input_dict::Dict{String, Union{Float64, Int}}: Input dictionary with parameter values
- required_params::Dict{String, Function}: Dictionary mapping parameter names to their 
    type conversion functions
- use_lists::Bool = true: If true, converts values to single-element vectors

# Returns
- Dict{String, Any}: Updated dictionary with converted values
"""
function update_initialization_dict(
    input_dict::Dict{String, Union{Float64, Int}},
    required_params::Dict{String, Function},
    use_lists::Bool = true
)::Dict{String, Any}

    updated_dict = Dict{String, Any}()

    if input_dict !== nothing
        required_names = keys(required_params)
        input_names = keys(input_dict)

        for required_name in required_names
            if !(required_name in input_names)
                error("initialization_params dict is missing parameter $required_name")
            end
            type_func = required_params[required_name]
            value = input_dict[required_name]
            if use_lists
                updated_dict[required_name] = [type_func(value)]
            else
                updated_dict[required_name] = type_func(value)
            end
        end
    end

    return updated_dict
end

"""
    check_dictionary_names(valid_keys::Vector{String}, 
        input_dict::Union{Dict{String, <:Union{Float64, Int64, String}}, Nothing})::Nothing

Validates that all keys in the input dictionary are present in the list of valid keys. 
Raises an error if any invalid keys are found.

# Arguments
- valid_keys::Vector{String}: List of valid key names
- input_dict::Union{Dict{String, <:Union{Float64, Int64, String}}, Nothing}: Input dictionary to validate

# Returns
- Nothing
"""
function check_dictionary_names(
    valid_keys::Vector{String},
    input_dict::Union{Dict{String, <:Union{Float64, Int64, String}}, Nothing}
)::Nothing
    if !isnothing(input_dict)
        for input_name in keys(input_dict)
            if !(input_name in valid_keys)
                error("Input key name $input_name is not a valid name for the input dictionary.")
            end
        end
    end
    return nothing
end

function check_for_missing_key_names(
    required_keys::Vector{String},
    input_dict::Dict{String, <:Union{Float64, Int, String}}
)::Nothing
    for key in required_keys
        if !haskey(input_dict, key)
            error("Missing key name: $key. Required keys are: $(required_keys)")
        end
    end
    return nothing
end

end # end module
