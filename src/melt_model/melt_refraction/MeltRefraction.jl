module MeltRefraction

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs

""" Update marker ID based on melt refraction.

# Updated Arrays
- model.markers.arrays.material.marker_matid.array
    - Material ID of each marker.
"""
function update_marker_type_for_melt_refraction!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    marker_matid = model.markers.arrays.material.marker_matid.array

    matid_types = model.materials.dicts.matid_types
    felsic_cc_fertile = matid_types["FelsicContinentalCrustFertile"]
    felsic_cc_pmolten = matid_types["FelsicContinentalCrustPartiallyMolten"]
    felsic_cc_refrac = matid_types["FelsicContinentalCrustRefactory"]
    solidified_granite = matid_types["SolidifiedGranite"]
    solidified_rhyolite = matid_types["SolidifiedRhyolite"]
    mafic_cc_fertile = matid_types["MaficContinentalCrustFertile"]
    mafic_cc_pmolten = matid_types["MaficContinentalCrustPartiallyMolten"]
    mafic_cc_refrac = matid_types["MaficContinentalCrustRefactory"]
    ultramafic_fertile = matid_types["UltramaficMantleFertile"]
    ultramafic_pmolten = matid_types["UltramaficMantlePartiallyMolten"]
    ultramafic_refrac = matid_types["UltramaficMantleRefactory"]

    marker_extractable_meltfrac = model.markers.arrays.melt.marker_extractable_meltfrac.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                extractable_melt_fraction = marker_extractable_meltfrac[imarker]
                material_id = marker_matid[imarker]
            end
            if extractable_melt_fraction < 0
                material_id = transform_marker_type_for_melt_refraction(
                    material_id,
                    felsic_cc_fertile,
                    felsic_cc_pmolten,
                    felsic_cc_refrac,
                    solidified_granite,
                    solidified_rhyolite,
                    mafic_cc_fertile,
                    mafic_cc_pmolten,
                    mafic_cc_refrac,
                    ultramafic_fertile,
                    ultramafic_pmolten,
                    ultramafic_refrac
                )
                @inbounds marker_matid[imarker] = material_id
            end
        end
    end
    return nothing
end

""" Transform material id based on melt refraction.

This is used to visualize zones of refractory materials.

# Arguments
- material_id::Int16
    - Current material ID of the marker.
- felsic_cc_fertile::Vector{Int16}
    - Material IDs for fertile felsic continental crust.
- felsic_cc_pmolten::Vector{Int16}
    - Material IDs for partially molten felsic continental crust.
- felsic_cc_refrac::Vector{Int16}
    - Material IDs for refractory felsic continental crust.
- solidified_granite::Vector{Int16}
    - Material IDs for solidified granite.
- solidified_rhyolite::Vector{Int16}
    - Material IDs for solidified rhyolite.
- mafic_cc_fertile::Vector{Int16}
    - Material IDs for fertile mafic continental crust.
- mafic_cc_pmolten::Vector{Int16}
    - Material IDs for partially molten mafic continental crust.
- mafic_cc_refrac::Vector{Int16}
    - Material IDs for refractory mafic continental crust.
- ultramafic_fertile::Vector{Int16}
    - Material IDs for fertile ultramafic mantle.
- ultramafic_pmolten::Vector{Int16}
    - Material IDs for partially molten ultramafic mantle.
- ultramafic_refrac::Vector{Int16}
    - Material IDs for refractory ultramafic mantle.

# Returns
- material_id_new::Int16
    - New material ID based on melt refraction.
"""
function transform_marker_type_for_melt_refraction(
    material_id::Int16,
    felsic_cc_fertile::Vector{Int16},
    felsic_cc_pmolten::Vector{Int16},
    felsic_cc_refrac::Vector{Int16},
    solidified_granite::Vector{Int16},
    solidified_rhyolite::Vector{Int16},
    mafic_cc_fertile::Vector{Int16},
    mafic_cc_pmolten::Vector{Int16},
    mafic_cc_refrac::Vector{Int16},
    ultramafic_fertile::Vector{Int16},
    ultramafic_pmolten::Vector{Int16},
    ultramafic_refrac::Vector{Int16}
)::Int16
    material_id_new = material_id
    if (
        material_id in felsic_cc_fertile || material_id in felsic_cc_pmolten ||
        material_id in solidified_granite || material_id in solidified_rhyolite
    )
        material_id_new = felsic_cc_refrac[1]
    elseif material_id in mafic_cc_fertile || material_id in mafic_cc_pmolten
        material_id_new = mafic_cc_refrac[1]
    elseif material_id in ultramafic_fertile || material_id in ultramafic_pmolten
        material_id_new = ultramafic_refrac[1]
    end
    return material_id_new
end

end # module 