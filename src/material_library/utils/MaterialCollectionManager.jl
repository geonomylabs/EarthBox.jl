module MaterialCollectionManager

"""
    MaterialCollection

# Fields
- `materials::NamedTuple`:
    - Named tuple of material names associated with the material collection.
- `path::String`:
    - Path to the yaml formatted material collection file with groups of materials and their properties.
"""
mutable struct MaterialCollection
    materials::NamedTuple
    path::String
end

end