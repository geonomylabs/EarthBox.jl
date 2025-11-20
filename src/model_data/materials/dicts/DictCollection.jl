module DictCollection

"""
    Dicts

Collection of material property dictionaries.

# Fields
- `matid_domains::Dict{String, Int64}`: Dictionary mapping domain name to a single material ID
- `matid_types::Dict{String, Vector{Int64}}`: Dictionary mapping type name to multiple material IDs

# Nested Dot Access
- `matid_domains = model.materials.dicts.matid_domains`
- `matid_types = model.materials.dicts.matid_types`

"""
mutable struct Dicts
    matid_domains::Dict{String, Int16}
    matid_types::Dict{String, Vector{Int16}}
end

end # module 