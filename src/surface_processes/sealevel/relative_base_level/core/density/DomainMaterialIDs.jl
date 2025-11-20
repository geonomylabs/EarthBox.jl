module DomainMaterialIDs

import EarthBox.ModelDataContainer: ModelData

function get_continental_crust_material_ids(model::ModelData)::Tuple{Int16,Int16}
    matid_domains = model.materials.dicts.matid_domains
    id_upper_continental_crust = matid_domains["UpperContinentalCrust"]
    id_lower_continental_crust = matid_domains["LowerContinentalCrust"]
    return id_upper_continental_crust, id_lower_continental_crust
end

function get_lithosphere_material_ids(model::ModelData)::Tuple{Int16,Int16,Int16}
    matid_domains = model.materials.dicts.matid_domains
    id_upper_mantle_lithosphere = matid_domains["UpperMantleLithosphere"]
    id_middle_mantle_lithosphere = matid_domains["MiddleMantleLithosphere"]
    id_lower_mantle_lithosphere = matid_domains["LowerMantleLithosphere"]
    return (
        id_upper_mantle_lithosphere,
        id_middle_mantle_lithosphere,
        id_lower_mantle_lithosphere
    )
end

function get_asthenosphere_material_id(model::ModelData)::Int16
    matid_domains = model.materials.dicts.matid_domains
    id_asthenosphere = matid_domains["Asthenosphere"]
    return id_asthenosphere
end

end # module