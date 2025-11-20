module GetExtractedMoltenIDs

import EarthBox.ModelDataContainer: ModelData

function get_extracted_gabbroic_molten_ids(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types

    if length(matid_types["SolidifiedGabbroPartiallyMolten"]) > 0
        matid_solidified_gabbro_partially_molten = matid_types[
        "SolidifiedGabbroPartiallyMolten"][1]
    else
        matid_solidified_gabbro_partially_molten = -1
    end
    if length(matid_types["ExtractedGabbroicMagma"]) > 0
        matid_extracted_gabbroic_magma = matid_types[
            "ExtractedGabbroicMagma"][1]
    else
        matid_extracted_gabbroic_magma = -1
    end
    if length(matid_types["SolidifiedLayeredGabbroPartiallyMolten"]) > 0
        matid_solidified_layered_gabbro_partially_molten = matid_types[
            "SolidifiedLayeredGabbroPartiallyMolten"][1]
    else
        matid_solidified_layered_gabbro_partially_molten = -1
    end
    if length(matid_types["ExtractedLayeredGabbroicMagma"]) > 0
        matid_extracted_layered_gabbroic_magma = matid_types[
            "ExtractedLayeredGabbroicMagma"][1]
    else
        matid_extracted_layered_gabbroic_magma = -1
    end

    return [
        matid_solidified_gabbro_partially_molten,
        matid_extracted_gabbroic_magma,
        matid_solidified_layered_gabbro_partially_molten,
        matid_extracted_layered_gabbroic_magma
    ]
end

function get_mantle_molten_ids(model::ModelData)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    matid_partially_molten_mantle = matid_types[
        "UltramaficMantlePartiallyMolten"][1]
    return [matid_partially_molten_mantle]
end

end # module 