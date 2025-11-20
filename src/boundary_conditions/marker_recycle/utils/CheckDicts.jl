module CheckDicts

"""
    check_domain_dict(domain_dict::Dict{String, Int}, dict_type::String, model_type::String)

Check if required material types or domains are present in the model.

# Arguments
- `domain_dict::Dict{String, Int}`: Dictionary of domain names or type names and associated material ids
- `dict_type::String`: Type of dictionary to check: either "type" or "domain"
- `model_type::String`: String describing model for which the dictionary is being checked

# Throws
- `ValueError`: If a required type or domain is not found in the materials input
"""
function check_domain_dict(
    domain_dict::Dict{String, Int16},
    dict_type::String,
    model_type::String
)::Nothing
    for (name, matid) in domain_dict
        if matid == -1
            throw(ErrorException(
                "$name $dict_type is required for $model_type " *
                "recycling model but not found in materials input."
            ))
        end
    end
    return nothing
end

end # module CheckDicts 