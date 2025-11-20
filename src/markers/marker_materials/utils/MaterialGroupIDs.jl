module MaterialGroupIDs

using EarthBox.ModelDataContainer: ModelData

"""
    get_ids_for_all_mantle_rocks(model::ModelData)::Vector{Int64}

Get material IDs for all mantle rocks.

This includes the following material types:
- UltramaficMantleFertile
- UltramaficMantlePartiallyMolten
- UltramaficMantleRefactory
"""
function get_ids_for_all_mantle_rocks(model::ModelData)::Vector{Int16}
    return get_all_mantle_ids_lists(model)
end

"""
    get_ids_for_asthenospheric_mantle(model::ModelData)::Vector{Int64}

Get material IDs for asthenospheric mantle.

This includes the following material types:
- UltramaficMantleFertile
- UltramaficMantlePartiallyMolten
- UltramaficMantleRefactory
"""
function get_ids_for_asthenospheric_mantle(model::ModelData)::Vector{Int16}
    material_id_lists = get_material_id_lists(model)
    return material_id_lists["asthenospheric_mantle"]
end

"""
    get_ids_for_felsic_and_mafic_continental_crust(model::ModelData)::Tuple{Vector{Int16}, Vector{Int16}}

Get material IDs for felsic and mafic continental crust.

Returns:
- felsic_continental_crust: Vector of felsic continental crust material ids
- mafic_continental_crust: Vector of mafic continental crust material ids
"""
function get_ids_for_felsic_and_mafic_continental_crust(
    model::ModelData
)::Tuple{Vector{Int16}, Vector{Int16}}
    material_id_lists = get_material_id_lists(model)
    return (material_id_lists["felsic_continental_crust"], material_id_lists["mafic_continental_crust"])
end

"""
    get_ids_for_plastic_materials_with_strain_weakening(model::ModelData)::Vector{Int16}

Get material IDs for frictional plastic materials with strain weakening.

This includes the following material types:
- weak_continental_crust: crust with strain weakening
- weak_mantle: mantle with strain weakening
- gabbroic_crust: gabbroic crust with strain weakening
- sediment: sediment with strain weakening
"""
function get_ids_for_plastic_materials_with_strain_weakening(model::ModelData)::Vector{Int16}
    material_id_lists = get_material_id_lists(model)
    return vcat(
        material_id_lists["weak_continental_crust"],
        material_id_lists["weak_mantle"],
        material_id_lists["gabbroic_crust"],
        material_id_lists["sediment"]
    )
end

"""
    get_ids_for_lithospheric_strong_zones(model::ModelData)::Vector{Int16}

Get material IDs for lithospheric strong zones.

This includes the following material types:
- strong_continental_crust: strong continental crust
- strong_mantle: strong mantle
"""
function get_ids_for_lithospheric_strong_zones(model::ModelData)::Vector{Int16}
    material_id_lists = get_material_id_lists(model)
    return vcat(
        material_id_lists["strong_continental_crust"],
        material_id_lists["strong_mantle"]
    )
end

function get_ids_for_solidified_basalt(model::ModelData)::Vector{Int16}
    material_id_lists = get_material_id_lists(model)
    return material_id_lists["solidified_basalt"]
end

"""
    get_ids_for_oceanic_crust(model::ModelData)::Vector{Int16}

Get material IDs for gabbroic basalt.

This includes solidified and partially molten gabbro and layered gabbro and
solidified basalt.
"""
function get_ids_for_oceanic_crust(model::ModelData)::Vector{Int16}
    material_id_lists = get_material_id_lists(model)
    return material_id_lists["gabbroic_crust"]
end

function get_ids_for_sediment(model::ModelData)::Vector{Int16}
    material_id_lists = get_material_id_lists(model)
    return material_id_lists["sediment"]
end

function get_serpentinite_ids_array(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    return collect(matid_types["Serpentinite"])
end

function get_mantle_ids_array(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    ids_mantle = collect(matid_types["UltramaficMantleFertile"])
    append!(ids_mantle, matid_types["UltramaficMantlePartiallyMolten"])
    append!(ids_mantle, matid_types["UltramaficMantleRefactory"])
    return ids_mantle
end

function get_gabbro_ids_array(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    ids_gabbro = collect(matid_types["SolidifiedGabbro"])
    append!(ids_gabbro, matid_types["SolidifiedLayeredGabbro"])
    return ids_gabbro
end

function get_magma_ids_array(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    ids_magma = collect(matid_types["ExtractedGabbroicMagma"])
    append!(ids_magma, matid_types["ExtractedLayeredGabbroicMagma"])
    return ids_magma
end

function get_continental_crust_ids_array(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    ids_continental_crust = collect(matid_types["FelsicContinentalCrustFertile"])
    append!(ids_continental_crust, matid_types["FelsicContinentalCrustPartiallyMolten"])
    append!(ids_continental_crust, matid_types["FelsicContinentalCrustRefactory"])
    append!(ids_continental_crust, matid_types["MaficContinentalCrustFertile"])
    append!(ids_continental_crust, matid_types["MaficContinentalCrustPartiallyMolten"])
    append!(ids_continental_crust, matid_types["MaficContinentalCrustRefactory"])
    return ids_continental_crust
end

function get_felsic_continental_crust_ids_array(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    ids_felsic_continental_crust = collect(matid_types["FelsicContinentalCrustFertile"])
    append!(ids_felsic_continental_crust, matid_types["FelsicContinentalCrustPartiallyMolten"])
    append!(ids_felsic_continental_crust, matid_types["FelsicContinentalCrustRefactory"])
    return ids_felsic_continental_crust
end

function get_mafic_continental_crust_ids_array(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    ids_mafic_continental_crust = collect(matid_types["MaficContinentalCrustFertile"])
    append!(ids_mafic_continental_crust, matid_types["MaficContinentalCrustPartiallyMolten"])
    append!(ids_mafic_continental_crust, matid_types["MaficContinentalCrustRefactory"])
    return ids_mafic_continental_crust
end

function get_strong_continental_crust_ids_array(model::ModelData)::Vector{Int16}
    matid_domains = model.materials.dicts.matid_domains
    id_upper_continental_crust_strong_zone = matid_domains["UpperContinentalCrustStrongZone"]
    id_lower_continental_crust_strong_zone = matid_domains["LowerContinentalCrustStrongZone"]
    return [id_upper_continental_crust_strong_zone, id_lower_continental_crust_strong_zone]
end

function get_mantle_strong_zone_id(model::ModelData)::Int16
    matid_domains = model.materials.dicts.matid_domains
    return matid_domains["LithosphericMantleStrongZone"]
end

function get_mantle_lithosphere_ids(model::ModelData)::Vector{Int16}
    matid_domains = model.materials.dicts.matid_domains
    id_upper_mantle_lithosphere = matid_domains["UpperMantleLithosphere"]
    id_middle_mantle_lithosphere = matid_domains["MiddleMantleLithosphere"]
    id_lower_mantle_lithosphere = matid_domains["LowerMantleLithosphere"]
    return [id_upper_mantle_lithosphere, id_middle_mantle_lithosphere, id_lower_mantle_lithosphere]
end

function get_solidified_basalt_ids(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    return [matid_types["SolidifiedBasalt"][1]]
end

""" Get dictionary of material id lists.

Returns a dictionary with the following keys:
- weak_continental_crust: list of normal continental crust material ids
- strong_continental_crust: list of strong continental crust material ids
- weak_mantle: list of normal mantle material ids
- strong_mantle: list of molten mantle material ids
- asthenospheric_mantle: list of asthenospheric mantle material ids
- felsic_continental_crust: list of felsic continental crust material ids
- mafic_continental_crust: list of mafic continental crust material ids
- gabbroic_crust: list of non-magma gabbroic crust material ids
- sediment: list of sediment material ids
- solidified_basalt: list of solidified basalt material ids
"""
function get_material_id_lists(model::ModelData)::Dict{String, Vector{Int16}}
    (
        ids_weak_continental_crust, ids_strong_continental_crust
    ) = get_continental_crust_id_lists(model)
    (
        ids_weak_mantle, ids_strong_mantle
    ) = get_mantle_ids_lists(model)
    ids_asthenospheric_mantle = get_asthenospheric_mantle_ids_lists(model)
    (
        ids_felsic_continental_crust, ids_mafic_continental_crust
    ) = get_continental_felsic_and_mafic_crust_ids_lists(model)
    ids_gabbroic_crust = get_non_magma_gabbroic_crust_id_list(model)
    ids_solidified_basalt = get_solidified_basalt_id_lists(model)

    return Dict{String, Vector{Int16}}(
        "weak_continental_crust" => ids_weak_continental_crust,
        "strong_continental_crust" => ids_strong_continental_crust,
        "weak_mantle" => ids_weak_mantle,
        "strong_mantle" => ids_strong_mantle,
        "asthenospheric_mantle" => ids_asthenospheric_mantle,
        "felsic_continental_crust" => ids_felsic_continental_crust,
        "mafic_continental_crust" => ids_mafic_continental_crust,
        "gabbroic_crust" => ids_gabbroic_crust,
        "sediment" => [get_sediment_id(model)],
        "solidified_basalt" => ids_solidified_basalt
    )
end

function get_solidified_basalt_id_lists(model::ModelData)::Vector{Int16}
    return collect(get_solidified_basalt_ids(model))
end

"""
    get_continental_crust_id_lists(model::ModelData)::Tuple{Vector{Int16}, Vector{Int16}}

Get array of continental crust material ids.

Returns:
- ids_weak_continental_crust: List of normal continental crust material ids
- ids_strong_continental_crust: List of strong continental crust material ids
"""
function get_continental_crust_id_lists(model::ModelData)::Tuple{Vector{Int16}, Vector{Int16}}
    strong_zone_ids = get_strong_continental_crust_ids_array(model)
    ids_continental_crust = get_continental_crust_ids_array(model)
    ids_weak_continental_crust = Int16[]
    ids_strong_continental_crust = Int16[]
    
    for matid in ids_continental_crust
        if matid ∉ strong_zone_ids
            push!(ids_weak_continental_crust, matid)
        else
            push!(ids_strong_continental_crust, matid)
        end
    end
    return (ids_weak_continental_crust, ids_strong_continental_crust)
end

"""
    get_mantle_ids_lists(model::ModelData)::Tuple{Vector{Int16}, Vector{Int16}}

Get array of mantle material ids.

Returns:
- ids_weak_mantle: List of normal mantle material ids
- ids_strong_mantle: List of molten mantle material ids
"""
function get_mantle_ids_lists(model::ModelData)::Tuple{Vector{Int16}, Vector{Int16}}
    ids_mantle = get_mantle_ids_array(model)
    id_mantle_strong_zone = get_mantle_strong_zone_id(model)
    ids_weak_mantle = Int16[]
    ids_strong_mantle = Int16[]
    
    for matid in ids_mantle
        if matid != id_mantle_strong_zone
            push!(ids_weak_mantle, matid)
        else
            push!(ids_strong_mantle, matid)
        end
    end
    return (ids_weak_mantle, ids_strong_mantle)
end

"""
    get_asthenospheric_mantle_ids_lists(model::ModelData)::Vector{Int16}

Get list of asthenospheric_mantle material ids.

Returns:
- ids_asthenospheric_mantle: List of asthenospheric mantle material ids
"""
function get_asthenospheric_mantle_ids_lists(model::ModelData)::Vector{Int16}
    ids_mantle = get_mantle_ids_array(model)
    ids_mantle_lithosphere = get_mantle_lithosphere_ids(model)
    ids_asthenospheric_mantle = Int16[]
    
    for matid in ids_mantle
        if matid ∉ ids_mantle_lithosphere
            push!(ids_asthenospheric_mantle, matid)
        end
    end
    return ids_asthenospheric_mantle
end

"""
    get_all_mantle_ids_lists(model::ModelData)::Vector{Int16}

Get list of asthenospheric_mantle material ids.

Returns:
- ids_mantle: List mantle material ids
"""
function get_all_mantle_ids_lists(model::ModelData)::Vector{Int16}
    return collect(get_mantle_ids_array(model))
end

"""
    get_continental_felsic_and_mafic_crust_ids_lists(model::ModelData)::Tuple{Vector{Int16}, Vector{Int16}}

Get lists of felsic and mafic continental crust material ids.

Returns:
- ids_felsic_continental_crust: List of felsic continental crust material ids
- ids_mafic_continental_crust: List of mafic continental crust material ids
"""
function get_continental_felsic_and_mafic_crust_ids_lists(
    model::ModelData
)::Tuple{Vector{Int16}, Vector{Int16}}
    ids_felsic_continental_crust = collect(get_felsic_continental_crust_ids_array(model))
    ids_mafic_continental_crust = collect(get_mafic_continental_crust_ids_array(model))
    return (ids_felsic_continental_crust, ids_mafic_continental_crust)
end

function get_non_magma_gabbroic_crust_id_list(model::ModelData)::Vector{Int16}
    return collect(get_non_magma_gabbroic_crust_ids_array(model))
end

function get_molten_gabbro_ids(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    matid_extracted_layered_gabbroic_magma = matid_types["ExtractedLayeredGabbroicMagma"][1]
    matid_solidified_layered_gabbro_partially_molten = matid_types["SolidifiedLayeredGabbroPartiallyMolten"][1]
    matid_extracted_gabbroic_magma = matid_types["ExtractedGabbroicMagma"][1]
    matid_solidified_gabbro_partially_molten = matid_types["SolidifiedGabbroPartiallyMolten"][1]
    
    return [
        matid_extracted_gabbroic_magma,
        matid_solidified_gabbro_partially_molten,
        matid_extracted_layered_gabbroic_magma,
        matid_solidified_layered_gabbro_partially_molten
    ]
end

function get_molten_layered_gabbro_ids(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    matid_extracted_layered_gabbroic_magma = matid_types["ExtractedLayeredGabbroicMagma"][1]
    matid_solidified_layered_gabbro_partially_molten = matid_types["SolidifiedLayeredGabbroPartiallyMolten"][1]
    
    return [
        matid_extracted_layered_gabbroic_magma,
        matid_solidified_layered_gabbro_partially_molten
    ]
end

function get_non_magma_gabbroic_crust_ids_array(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    matid_solidified_gabbro = matid_types["SolidifiedGabbro"][1]
    matid_solidified_gabbro_partially_molten = matid_types["SolidifiedGabbroPartiallyMolten"][1]
    matid_solidified_layered_gabbro = matid_types["SolidifiedLayeredGabbro"][1]
    matid_solidified_layered_gabbro_partially_molten = matid_types["SolidifiedLayeredGabbroPartiallyMolten"][1]
    matid_solidified_basalt = matid_types["SolidifiedBasalt"][1]
    
    return [
        matid_solidified_gabbro,
        matid_solidified_gabbro_partially_molten,
        matid_solidified_layered_gabbro,
        matid_solidified_layered_gabbro_partially_molten,
        matid_solidified_basalt
    ]
end

function get_sediment_id(model::ModelData)::Int16
    matid_types = model.materials.dicts.matid_types
    matid_array = matid_types["Sediment"]
    if isempty(matid_array)
        return -1
    else
        return matid_array[1]
    end
end

function get_sticky_material_ids(model::ModelData)::Tuple{Int16, Int16}
    types = model.materials.dicts.matid_types
    matid_sticky_air = types["StickyAir"][1]
    matid_sticky_water = types["StickyWater"][1]
    return (matid_sticky_air, matid_sticky_water)
end

end # module MaterialGroupIDs 