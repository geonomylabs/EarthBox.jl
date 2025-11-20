module LoadMaterialDicts

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import ..MaterialContainer: Material
import ..MaterialContainer: get_material_domains, MaterialDomainsRegistry
import ..MaterialContainer: get_material_types
import ..LoadMaterialArrays: MaterialCollectionDictType

""" Update initialized material ID dictionaries using material model.
"""
function update_material_id_dicts(
    materials::MaterialCollectionDictType, 
    model::ModelData
)::Nothing
    domain_matid_dict, type_matid_dict = get_material_id_dicts(materials)
    for (key, value) in domain_matid_dict
        model.materials.dicts.matid_domains[key] = value
    end
    for (key, value) in type_matid_dict
        model.materials.dicts.matid_types[key] = value
    end
end

"""Get material id dictionaries for material domain and material type.
Each domain is associated with a single material id.
"""
function get_material_id_dicts(
    materials::MaterialCollectionDictType
)::Tuple{Dict{String, Int64}, Dict{String, Vector{Int16}}}
    domain_matid_dict = define_material_domain_ids(materials)
    type_matid_dict = define_material_type_ids(materials)
    return domain_matid_dict, type_matid_dict
end

""" Define material id for each domain.
"""
function define_material_domain_ids(
    materials::MaterialCollectionDictType
)::Dict{String, Int64}
    domain_matid_dict = Dict{String, Int64}()
    mat_domains = get_material_domains()
    for mat_domain in mat_domains
        domain_matid_dict[mat_domain] = get_matid_for_domain(materials, mat_domain)
    end

    mdr = MaterialDomainsRegistry()
    if domain_matid_dict[mdr.mantle_lithosphere] >= 1
        mld = domain_matid_dict[mdr.mantle_lithosphere]
        domain_matid_dict[mdr.upper_mantle_lithosphere] = mld
        domain_matid_dict[mdr.middle_mantle_lithosphere] = mld
        domain_matid_dict[mdr.lower_mantle_lithosphere] = mld
        print_info("Mantle lithosphere domain $(mld) has been assigned to upper, middle, and lower mantle lithosphere domains.", level=1)
    end

    return domain_matid_dict
end

""" Get unique integer material id for an input material domain string.

Loop over all material objects associated with model and assign a
material ID if material objects domain is equal to the input domain.

If multiple material objects have the same domain, then the last
object in the materials dictionary will be used to assign the ID. That
is each domain can have only one ID.

# Arguments
- `materials`: Dictionary defining the relationship between material ID's and material objects
- `mat_domain`: Name of material domain for which a material id needs to be assigned

# Returns
- `matid_domain`: Material id associated with material domain
"""
function get_matid_for_domain(
    materials::MaterialCollectionDictType, 
    mat_domain::String
)::Int64
    matid_domain = -1
    for (matid, mat_obj) in materials
        mat_domain_check = mat_obj.mat_domain.value
        if mat_domain_check == mat_domain
            if !isa(matid, Int16)
                matid_domain = parse(Int16, matid)
            else
                matid_domain = matid
            end
        end
    end
    return matid_domain
end

""" Define material id's associated with each type.
"""
function define_material_type_ids(materials::MaterialCollectionDictType)::Dict{String, Vector{Int16}}
    type_matid_dict = Dict{String, Vector{Int16}}()
    mat_types = get_material_types()
    for mat_type in mat_types
        matid_array = get_matid_array_for_type(materials, mat_type)
        type_matid_dict[mat_type] = matid_array
    end
    return type_matid_dict
end

""" Get array of material id's associated with material type.
"""
function get_matid_array_for_type(
    materials::MaterialCollectionDictType, 
    mat_type::String
)::Vector{Int16}
    matid_type_list = Int16[]
    for (matid, mat_obj) in materials
        mat_type_check = mat_obj.mat_type.value
        if mat_type_check == mat_type
            if isa(matid, String)
                push!(matid_type_list, parse(Int16, matid))
            else
                push!(matid_type_list, matid)
            end
        end
    end
    return matid_type_list
end

end # module