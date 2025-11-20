module GetMaterialIDs

import ..MaterialsContainer: Materials
import ..MaterialsContainer.LoadMaterialDicts: get_material_id_dicts
import ..MaterialsContainer.LoadMaterialDicts: define_material_domain_ids
import ..MaterialsContainer.LoadMaterialDicts: define_material_type_ids

""" Get material IDs for continental crust.
"""
function get_ids_for_continental_crust(materials::Materials)::Vector{Int16}
    (
        material_ids_felsic_continental_crust,
        material_ids_mafic_continental_crust
    ) = get_ids_for_felsic_and_mafic_continental_crust(materials)
    material_ids_continental_crust = vcat(
        material_ids_felsic_continental_crust,
        material_ids_mafic_continental_crust
    )
    return material_ids_continental_crust
end

""" Get material IDs for asthenospheric mantle.

This includes the following material types:
- UltramaficMantleFertile
- UltramaficMantlePartiallyMolten
- UltramaficMantleRefactory
"""
function get_ids_for_asthenospheric_mantle(materials::Materials)::Vector{Int16}
    material_id_lists = get_material_id_lists(materials)
    asthenospheric_mantle = material_id_lists["asthenospheric_mantle"]
    return asthenospheric_mantle
end

""" Get material IDs for felsic and mafic continental crust.

Returns
-------
felsic_continental_crust : Vector{Int16}
    List of felsic continental crust material ids.

mafic_continental_crust : Vector{Int16}
    List of mafic continental crust material ids.
"""
function get_ids_for_felsic_and_mafic_continental_crust(
    materials::Materials
)::Tuple{Vector{Int16}, Vector{Int16}}
    material_id_lists = get_material_id_lists(materials)
    felsic_continental_crust = material_id_lists["felsic_continental_crust"]
    mafic_continental_crust = material_id_lists["mafic_continental_crust"]
    return felsic_continental_crust, mafic_continental_crust
end

""" Get material IDs for frictional plastic materials with strain weakening.

This includes the following material types:
- weak_continental_crust: crust with strain weakening
- weak_mantle: mantle with strain weakening
- gabbroic_crust: gabbroic crust with strain weakening
- sediment: sediment with strain weakening
"""
function get_ids_for_plastic_materials_with_strain_weakening(
    materials::Materials
)::Vector{Int16}
    material_id_lists = get_material_id_lists(materials)

    weak_continental_crust = material_id_lists["weak_continental_crust"]
    weak_mantle = material_id_lists["weak_mantle"]
    gabbroic_crust = material_id_lists["gabbroic_crust"]
    sediment = material_id_lists["sediment"]

    material_ids_to_update_plasticity = vcat(
        weak_continental_crust,
        weak_mantle,
        gabbroic_crust,
        sediment
    )

    return material_ids_to_update_plasticity
end

""" Get material IDs for lithospheric strong zones.

This includes the following material types:
- strong_continental_crust: strong continental crust
- strong_mantle: strong mantle
"""
function get_ids_for_lithospheric_strong_zones(
    materials::Materials
)::Vector{Int16}
    material_id_lists = get_material_id_lists(materials)

    strong_continental_crust = material_id_lists["strong_continental_crust"]
    strong_mantle = material_id_lists["strong_mantle"]

    material_ids_strong_zone = vcat(strong_continental_crust, strong_mantle)

    return material_ids_strong_zone
end

""" Get array of serpentinite material ids.
"""
function get_serpentinite_ids_array(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)
    ids_serpentinite = collect(matid_types["Serpentinite"])
    return ids_serpentinite
end

""" Get array of mantle material ids.
"""
function get_mantle_ids_array(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)
    ids_mantle = collect(matid_types["UltramaficMantleFertile"])
    append!(ids_mantle, matid_types["UltramaficMantlePartiallyMolten"])
    append!(ids_mantle, matid_types["UltramaficMantleRefactory"])
    return ids_mantle
end

""" Get array of gabbro material ids.
"""
function get_gabbro_ids_array(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)
    ids_gabbro = collect(matid_types["SolidifiedGabbro"])
    append!(ids_gabbro, matid_types["SolidifiedLayeredGabbro"])
    return ids_gabbro
end

""" Get array of magma material ids.
"""
function get_magma_ids_array(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)
    ids_magma = collect(matid_types["ExtractedGabbroicMagma"])
    append!(ids_magma, matid_types["ExtractedLayeredGabbroicMagma"])
    return ids_magma
end

""" Get array of continental crust material ids.
"""
function get_continental_crust_ids_array(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)
    ids_continental_crust = collect(matid_types["FelsicContinentalCrustFertile"])
    append!(ids_continental_crust, matid_types["FelsicContinentalCrustPartiallyMolten"])
    append!(ids_continental_crust, matid_types["FelsicContinentalCrustRefactory"])
    append!(ids_continental_crust, matid_types["MaficContinentalCrustFertile"])
    append!(ids_continental_crust, matid_types["MaficContinentalCrustPartiallyMolten"])
    append!(ids_continental_crust, matid_types["MaficContinentalCrustRefactory"])
    return ids_continental_crust
end

""" Get array of felsic continental crust material ids.
"""
function get_felsic_continental_crust_ids_array(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)
    ids_felsic_continental_crust = collect(matid_types["FelsicContinentalCrustFertile"])
    append!(ids_felsic_continental_crust, matid_types["FelsicContinentalCrustPartiallyMolten"])
    append!(ids_felsic_continental_crust, matid_types["FelsicContinentalCrustRefactory"])
    return ids_felsic_continental_crust
end

""" Get array of mafic continental crust material ids.
"""
function get_mafic_continental_crust_ids_array(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)
    ids_mafic_continental_crust = collect(matid_types["MaficContinentalCrustFertile"])
    append!(ids_mafic_continental_crust, matid_types["MaficContinentalCrustPartiallyMolten"])
    append!(ids_mafic_continental_crust, matid_types["MaficContinentalCrustRefactory"])
    return ids_mafic_continental_crust
end

""" Get array of strong continental crust material ids.
"""
function get_strong_continental_crust_ids_array(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    matid_domains, _matid_types = get_material_id_dicts(materials_dict)
    id_upper_continental_crust_strong_zone = matid_domains["UpperContinentalCrustStrongZone"]
    id_lower_continental_crust_strong_zone = matid_domains["LowerContinentalCrustStrongZone"]
    return [id_upper_continental_crust_strong_zone, id_lower_continental_crust_strong_zone]
end

""" Get mantle lithosphere strong zone material ids.
"""
function get_mantle_strong_zone_id(materials::Materials)::Int16
    materials_dict = materials.materials
    matid_domains, _matid_types = get_material_id_dicts(materials_dict)
    id_upper_mantle_lithosphere_strong_zone = matid_domains["LithosphericMantleStrongZone"]
    return id_upper_mantle_lithosphere_strong_zone
end

""" Get mantle lithosphere material ids.
"""
function get_mantle_lithosphere_ids(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    matid_domains, _matid_types = get_material_id_dicts(materials_dict)
    id_upper_mantle_lithosphere = matid_domains["UpperMantleLithosphere"]
    id_middle_mantle_lithosphere = matid_domains["MiddleMantleLithosphere"]
    id_lower_mantle_lithosphere = matid_domains["LowerMantleLithosphere"]
    return [id_upper_mantle_lithosphere, id_middle_mantle_lithosphere, id_lower_mantle_lithosphere]
end

""" Get dictionary of material id lists.

Returns
-------
material_id_lists : Dict{String, Vector{Int16}}
    Dictionary of material id lists:
    - 'weak_continental_crust': list of normal continental crust material ids.
    - 'strong_continental_crust': list of strong continental crust material ids.
    - 'weak_mantle': list of normal mantle material ids.
    - 'strong_mantle': list of molten mantle material ids.
    - 'asthenospheric_mantle': list of asthenospheric mantle material ids.
    - 'felsic_continental_crust': list of felsic continental crust material ids.
    - 'mafic_continental_crust': list of mafic continental crust material ids.
    - 'gabbroic_crust': list of non-magma gabbroic crust material ids.
    - 'sediment': list of sediment material ids.
"""
function get_material_id_lists(materials::Materials)::Dict{String, Vector{Int16}}
    (
        ids_weak_continental_crust,
        ids_strong_continental_crust
    ) = get_continental_crust_id_lists(materials)

    (
        ids_weak_mantle,
        ids_strong_mantle
    ) = get_mantle_ids_lists(materials)

    ids_asthenospheric_mantle = get_asthenospheric_mantle_ids_lists(materials)

    (
        ids_felsic_continental_crust,
        ids_mafic_continental_crust
    ) = get_continental_felsic_and_mafic_crust_ids_lists(materials)

    ids_gabbroic_crust = get_non_magma_gabbroic_crust_id_list(materials)

    material_id_lists = Dict{String, Vector{Int16}}(
        "weak_continental_crust" => ids_weak_continental_crust,
        "strong_continental_crust" => ids_strong_continental_crust,
        "weak_mantle" => ids_weak_mantle,
        "strong_mantle" => ids_strong_mantle,
        "asthenospheric_mantle" => ids_asthenospheric_mantle,
        "felsic_continental_crust" => ids_felsic_continental_crust,
        "mafic_continental_crust" => ids_mafic_continental_crust,
        "gabbroic_crust" => ids_gabbroic_crust,
        "sediment" => [get_sediment_id(materials)]
    )

    return material_id_lists
end

""" Get array of continental crust material ids.

Returns
-------
ids_weak_continental_crust : Vector{Int16}
    List of normal continental crust material ids.

ids_strong_continental_crust : Vector{Int16}
    List of strong continental crust material ids.
"""
function get_continental_crust_id_lists(
    materials::Materials
)::Tuple{Vector{Int16}, Vector{Int16}}
    strong_zone_ids = get_strong_continental_crust_ids_array(materials)
    ids_continental_crust = get_continental_crust_ids_array(materials)
    ids_weak_continental_crust = Int64[]
    ids_strong_continental_crust = Int64[]
    for matid in ids_continental_crust
        if matid ∉ strong_zone_ids
            push!(ids_weak_continental_crust, matid)
        else
            push!(ids_strong_continental_crust, matid)
        end
    end
    return ids_weak_continental_crust, ids_strong_continental_crust
end

""" Get array of mantle material ids.

Returns
-------
ids_weak_mantle : Vector{Int16}
    List of normal mantle material ids.

ids_strong_mantle : Vector{Int16}
    List of molten mantle material ids.
"""
function get_mantle_ids_lists(
    materials::Materials
)::Tuple{Vector{Int16}, Vector{Int16}}
    ids_mantle = get_mantle_ids_array(materials)
    id_mantle_strong_zone = get_mantle_strong_zone_id(materials)
    ids_weak_mantle = Int64[]
    ids_strong_mantle = Int64[]
    for matid in ids_mantle
        if matid != id_mantle_strong_zone
            push!(ids_weak_mantle, matid)
        else
            push!(ids_strong_mantle, matid)
        end
    end
    return ids_weak_mantle, ids_strong_mantle
end

""" Get list of asthenospheric_mantle material ids.

Returns
-------
ids_asthenospheric_mantle : Vector{Int16}
    List of asthenospheric mantle material ids.
"""
function get_asthenospheric_mantle_ids_lists(
    materials::Materials
)::Vector{Int16}
    ids_mantle = get_mantle_ids_array(materials)
    ids_mantle_lithosphere = get_mantle_lithosphere_ids(materials)

    ids_asthenospheric_mantle = Int64[]
    for matid in ids_mantle
        if matid ∉ ids_mantle_lithosphere
            push!(ids_asthenospheric_mantle, matid)
        end
    end
    return ids_asthenospheric_mantle
end

""" Get lists of felsic and mafic continental crust material ids.

Returns
-------
ids_felsic_continental_crust : Vector{Int16}
    List of felsic continental crust material ids.

ids_mafic_continental_crust : Vector{Int16}
    List of mafic continental crust material ids.
"""
function get_continental_felsic_and_mafic_crust_ids_lists(
    materials::Materials
)::Tuple{Vector{Int16}, Vector{Int16}}
    ids_felsic_continental_crust = get_felsic_continental_crust_ids_array(materials)
    ids_mafic_continental_crust = get_mafic_continental_crust_ids_array(materials)
    return collect(ids_felsic_continental_crust), collect(ids_mafic_continental_crust)
end

""" Get list of non-magma gabbroic crust material ids.
"""
function get_non_magma_gabbroic_crust_id_list(
    materials::Materials
)::Vector{Int16}
    ids_gabbro = get_non_magma_gabbroic_crust_ids_array(materials)
    return collect(ids_gabbro)
end

""" Get the molten gabbro material IDs.
"""
function get_molten_gabbro_ids(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)

    matid_extracted_layered_gabbroic_magma = matid_types["ExtractedLayeredGabbroicMagma"][1]
    matid_solidified_layered_gabbro_partially_molten = matid_types["SolidifiedLayeredGabbroPartiallyMolten"][1]
    matid_extracted_gabbroic_magma = matid_types["ExtractedGabbroicMagma"][1]
    matid_solidified_gabbro_partially_molten = matid_types["SolidifiedGabbroPartiallyMolten"][1]

    molten_gabbro_ids = [
        matid_extracted_gabbroic_magma,
        matid_solidified_gabbro_partially_molten,
        matid_extracted_layered_gabbroic_magma,
        matid_solidified_layered_gabbro_partially_molten
    ]
    return molten_gabbro_ids
end

""" Get the molten layered gabbro material IDs.
"""
function get_molten_layered_gabbro_ids(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)
    matid_extracted_layered_gabbroic_magma = matid_types["ExtractedLayeredGabbroicMagma"][1]
    matid_solidified_layered_gabbro_partially_molten = matid_types["SolidifiedLayeredGabbroPartiallyMolten"][1]

    molten_gabbro_ids = [
        matid_extracted_layered_gabbroic_magma,
        matid_solidified_layered_gabbro_partially_molten
    ]
    return molten_gabbro_ids
end

""" Get the molten gabbro material IDs.
"""
function get_non_magma_gabbroic_crust_ids_array(materials::Materials)::Vector{Int16}
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)
    matid_solidified_gabbro = matid_types["SolidifiedGabbro"][1]
    matid_solidified_gabbro_partially_molten = matid_types["SolidifiedGabbroPartiallyMolten"][1]
    matid_solidified_layered_gabbro = matid_types["SolidifiedLayeredGabbro"][1]
    matid_solidified_layered_gabbro_partially_molten = matid_types["SolidifiedLayeredGabbroPartiallyMolten"][1]
    matid_solidified_basalt = matid_types["SolidifiedBasalt"][1]

    non_magma_gabbro_ids = [
        matid_solidified_gabbro,
        matid_solidified_gabbro_partially_molten,
        matid_solidified_layered_gabbro,
        matid_solidified_layered_gabbro_partially_molten,
        matid_solidified_basalt
    ]
    return non_magma_gabbro_ids
end

""" Get sticky air/water material IDs tuple.
"""
function get_sticky_material_ids(materials::Materials)::Tuple{Int16, Int16}
    types = get_material_type_dict(materials)
    id_array = types["StickyAir"]
    if isempty(id_array)
        matid_sticky_air = Int16(-1)
    else
        matid_sticky_air = id_array[1]
    end
    id_array = types["StickyWater"]
    if isempty(id_array)
        matid_sticky_water = Int16(-1)
    else
        matid_sticky_water = id_array[1]
    end
    return (matid_sticky_air, matid_sticky_water)
end

""" Get the sediment material ID.
"""
function get_sediment_id(materials::Materials)::Int16
    materials_dict = materials.materials
    _matid_domains, matid_types = get_material_id_dicts(materials_dict)
    matid_array = matid_types["Sediment"]
    if isempty(matid_array)
        matid_sediment = Int16(-1)
    else
        matid_sediment = matid_array[1]
    end
    return matid_sediment
end

""" Get sediment material id.
"""
function get_sediment_material_id(materials::Materials)::Int16
    types = get_material_type_dict(materials)
    id_array = types["Sediment"]
    if !isempty(id_array)
        matid_sediment = types["Sediment"][1]
    else
        matid_sediment = Int16(-1)
    end
    return matid_sediment
end

""" Get solidified basalt material id.
"""
function get_solidified_basalt_material_id(materials::Materials)::Int16
    types = get_material_type_dict(materials)
    id_array = types["SolidifiedBasalt"]
    if !isempty(id_array)
        matid_basalt = types["SolidifiedBasalt"][1]
    else
        matid_basalt = Int16(-1)
    end
    return matid_basalt
end

""" Get layered gabbro material id.
"""
function get_solidified_layered_gabbro_material_id(materials::Materials)::Int16
    types = get_material_type_dict(materials)
    id_array = types["SolidifiedLayeredGabbro"]
    if !isempty(id_array)
        matid_layered_gabbro = types["SolidifiedLayeredGabbro"][1]
    else
        matid_layered_gabbro = Int16(-1)
    end
    return matid_layered_gabbro
end

""" Get gabbro material id.
"""
function get_solidified_gabbro_material_id(materials::Materials)::Int16
    types = get_material_type_dict(materials)
    id_array = types["SolidifiedGabbro"]
    if !isempty(id_array)
        matid_gabbro = types["SolidifiedGabbro"][1]
    else
        matid_gabbro = Int16(-1)
    end
    return matid_gabbro
end

""" Get salt material id.
"""
function get_salt_material_id(materials::Materials)::Int16
    types = get_material_type_dict(materials)
    matid_salt = types["Salt"][1]
    return matid_salt
end

""" Get mantle melting material ids.
"""
function get_mantle_melting_matids(materials::Materials)::Vector{Int16}
    types = get_material_type_dict(materials)

    fertile = types["UltramaficMantleFertile"]
    nfertile = length(fertile)

    partially_molten = types["UltramaficMantlePartiallyMolten"]
    npartially_molten = length(partially_molten)

    refractory = types["UltramaficMantleRefactory"]
    nrefractory = length(refractory)

    ntotal = nfertile + npartially_molten + nrefractory

    mantle_melting_mat_ids = zeros(Int16, ntotal)
    icount = 1
    for i in 1:nfertile
        mantle_melting_mat_ids[icount] = fertile[i]
        icount += 1
    end
    for i in 1:npartially_molten
        mantle_melting_mat_ids[icount] = partially_molten[i]
        icount += 1
    end
    for i in 1:nrefractory
        mantle_melting_mat_ids[icount] = refractory[i]
        icount += 1
    end
    return mantle_melting_mat_ids
end

""" Get mantle and gabbroic melting materials id's.
"""
function get_mantle_and_gabbroic_melting_matids(materials::Materials)::Vector{Int16}
    types = get_material_type_dict(materials)

    fertile = types["UltramaficMantleFertile"]
    nfertile = length(fertile)

    partially_molten = types["UltramaficMantlePartiallyMolten"]
    npartially_molten = length(partially_molten)

    refractory = types["UltramaficMantleRefactory"]
    nrefractory = length(refractory)

    gabbro_partially_molten = types["SolidifiedGabbroPartiallyMolten"]
    ngabbro_partially_molten = length(gabbro_partially_molten)

    layered_gabbro_partially_molten = types["SolidifiedLayeredGabbroPartiallyMolten"]
    nlayered_gabbro_partially_molten = length(layered_gabbro_partially_molten)

    ntotal = (
        nfertile +
        npartially_molten +
        nrefractory +
        ngabbro_partially_molten +
        nlayered_gabbro_partially_molten
    )

    mantle_and_gabbroic_melting_mat_ids = zeros(Int16, ntotal)
    icount = 1
    for i in 1:nfertile
        mantle_and_gabbroic_melting_mat_ids[icount] = fertile[i]
        icount += 1
    end
    for i in 1:npartially_molten
        mantle_and_gabbroic_melting_mat_ids[icount] = partially_molten[i]
        icount += 1
    end
    for i in 1:nrefractory
        mantle_and_gabbroic_melting_mat_ids[icount] = refractory[i]
        icount += 1
    end
    for i in 1:ngabbro_partially_molten
        mantle_and_gabbroic_melting_mat_ids[icount] = gabbro_partially_molten[i]
        icount += 1
    end
    for i in 1:nlayered_gabbro_partially_molten
        mantle_and_gabbroic_melting_mat_ids[icount] = layered_gabbro_partially_molten[i]
        icount += 1
    end
    return mantle_and_gabbroic_melting_mat_ids
end

function get_material_domain_dict(materials::Materials)::Dict{String, Int16}
    domain_matid_dict = define_material_domain_ids(materials.materials)
    return domain_matid_dict
end

function get_material_type_dict(materials::Materials)::Dict{String, Vector{Int16}}
    type_matid_dict = define_material_type_ids(materials.materials)
    return type_matid_dict
end

end 