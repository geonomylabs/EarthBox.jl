module DictCollection

"""
    Dicts

Collection of material property dictionaries.

# Fields
- `matid_domains::Dict{String, Int64}`: Dictionary mapping domain name to a single material ID
- `matid_types::Dict{String, Vector{Int64}}`: Dictionary mapping type name to multiple material IDs
- `cached_oceanic_crust_ids::Vector{Int16}`: Lazy-cached concatenation of
    oceanic-crust matids used by `Fractionation.get_matids_oceanic_crust`.
    Empty until first access; populated from `matid_types` on first call
    after materials initialization.

# Nested Dot Access
- `matid_domains = model.materials.dicts.matid_domains`
- `matid_types = model.materials.dicts.matid_types`
- `cached_oceanic_crust_ids = model.materials.dicts.cached_oceanic_crust_ids`

"""
mutable struct Dicts
    matid_domains::Dict{String, Int16}
    matid_types::Dict{String, Vector{Int16}}
    cached_oceanic_crust_ids::Vector{Int16}
end

end # module 