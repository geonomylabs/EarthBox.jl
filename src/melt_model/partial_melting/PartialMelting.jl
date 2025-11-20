module PartialMelting

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs

""" Update marker ID and type based on partial melting.

This function loops over markers and calls functions that update the marker
material id array ``marker_matid`` to account for partial melting of
solidified material when melt fraction is between 0 and 1.

# Updated Arrays
- `model.markers.arrays.material.marker_matid`: Material ID of the marker, 
  updated to reflect partial melting.
"""
function update_solidified_marker_type_for_partial_melting!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array

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

    solidified_gabbro = matid_types["SolidifiedGabbro"]
    solidified_gabbro_pmolten = matid_types["SolidifiedGabbroPartiallyMolten"]

    solidified_layered_gabbro = matid_types["SolidifiedLayeredGabbro"]
    solidified_layered_gabbro_pmolten = matid_types["SolidifiedLayeredGabbroPartiallyMolten"]

    solidified_basalt = matid_types["SolidifiedBasalt"]

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            melt_fraction = marker_meltfrac[imarker]
            material_id = marker_matid[imarker]
            if 0 < melt_fraction < 1
                material_id = transform_solid_marker_to_partially_molten(
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
                    ultramafic_refrac,
                    solidified_gabbro,
                    solidified_gabbro_pmolten,
                    solidified_basalt,
                    solidified_layered_gabbro,
                    solidified_layered_gabbro_pmolten
                )
            elseif melt_fraction <= 0
                material_id = transform_partially_molten_marker_to_zero_melt_state(
                    material_id,
                    felsic_cc_fertile,
                    felsic_cc_pmolten,
                    mafic_cc_fertile,
                    mafic_cc_pmolten,
                    ultramafic_fertile,
                    ultramafic_pmolten,
                    solidified_gabbro,
                    solidified_gabbro_pmolten,
                    solidified_layered_gabbro,
                    solidified_layered_gabbro_pmolten
                )
            end
            marker_matid[imarker] = material_id
        end
    end
    return nothing
end

""" Check if marker is magma.

# Arguments
- `material_id::Int16`: Material ID of the marker
- `extracted_gabbroic_magma::Vector{Int16}`: Array of extracted gabbroic magma IDs
- `extracted_layered_gabbroic_magma::Vector{Int16}`: Array of extracted layered 
  gabbroic magma IDs
- `extruded_gabbroic_magma::Vector{Int16}`: Array of extruded gabbroic magma IDs

# Returns
- `Bool`: True if marker is magma, False otherwise
"""
@inline function check_for_magma_marker(
    material_id::Int16,
    extracted_gabbroic_magma::Vector{Int16},
    extracted_layered_gabbroic_magma::Vector{Int16},
    extruded_gabbroic_magma::Vector{Int16}
)::Bool
    return material_id in extracted_gabbroic_magma ||
           material_id in extracted_layered_gabbroic_magma ||
           material_id in extruded_gabbroic_magma
end

""" Transform material id based on melt fraction.

This is used to visualize zones of partial melt.

# Arguments
- `material_id::Int16`: Material ID of the marker
- `extracted_gabbroic_magma::Vector{Int16}`: Array of extracted gabbroic magma IDs
- `extracted_layered_gabbroic_magma::Vector{Int16}`: Array of extracted layered 
  gabbroic magma IDs
- `extruded_gabbroic_magma::Vector{Int16}`: Array of extruded gabbroic magma IDs
- `solidified_gabbro_pmolten::Vector{Int16}`: Array of solidified gabbro 
  partially molten IDs
- `solidified_layered_gabbro_pmolten::Vector{Int16}`: Array of solidified layered 
  gabbro partially molten IDs

# Returns
- `Int16`: New material ID
"""
@inline function transform_magma_marker_to_partially_molten(
    material_id::Int16,
    extracted_gabbroic_magma::Vector{Int16},
    extracted_layered_gabbroic_magma::Vector{Int16},
    extruded_gabbroic_magma::Vector{Int16},
    solidified_gabbro_pmolten::Vector{Int16},
    solidified_layered_gabbro_pmolten::Vector{Int16}
)::Int16
    material_id_new = material_id
    if material_id in extracted_gabbroic_magma
        material_id_new = solidified_gabbro_pmolten[1]
    elseif material_id in extracted_layered_gabbroic_magma
        material_id_new = solidified_layered_gabbro_pmolten[1]
    elseif material_id in extruded_gabbroic_magma
        material_id_new = solidified_gabbro_pmolten[1]
    end
    return material_id_new
end

""" Transform material id based on melt fraction.

This is used to visualize zones of partial melt.

# Arguments
- `material_id::Int16`: Material ID of the marker
- `felsic_cc_fertile::Vector{Int16}`: Array of felsic continental crust fertile IDs
- `felsic_cc_pmolten::Vector{Int16}`: Array of felsic continental crust partially 
  molten IDs
- `felsic_cc_refrac::Vector{Int16}`: Array of felsic continental crust refractory IDs
- `solidified_granite::Vector{Int16}`: Array of solidified granite IDs
- `solidified_rhyolite::Vector{Int16}`: Array of solidified rhyolite IDs
- `mafic_cc_fertile::Vector{Int16}`: Array of mafic continental crust fertile IDs
- `mafic_cc_pmolten::Vector{Int16}`: Array of mafic continental crust partially 
  molten IDs
- `mafic_cc_refrac::Vector{Int16}`: Array of mafic continental crust refractory IDs
- `ultramafic_fertile::Vector{Int16}`: Array of ultramafic mantle fertile IDs
- `ultramafic_pmolten::Vector{Int16}`: Array of ultramafic mantle partially 
  molten IDs
- `ultramafic_refrac::Vector{Int16}`: Array of ultramafic mantle refractory IDs
- `solidified_gabbro::Vector{Int16}`: Array of solidified gabbro IDs
- `gabbro_pmolten::Vector{Int16}`: Array of gabbro partially molten IDs
- `solidified_basalt::Vector{Int16}`: Array of solidified basalt IDs
- `solidified_layered_gabbro::Vector{Int16}`: Array of solidified layered gabbro IDs
- `solidified_layered_gabbro_pmolten::Vector{Int16}`: Array of solidified layered 
  gabbro partially molten IDs

# Returns
- `Int16`: New material ID
"""
@inline function transform_solid_marker_to_partially_molten(
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
    ultramafic_refrac::Vector{Int16},
    solidified_gabbro::Vector{Int16},
    gabbro_pmolten::Vector{Int16},
    solidified_basalt::Vector{Int16},
    solidified_layered_gabbro::Vector{Int16},
    solidified_layered_gabbro_pmolten::Vector{Int16}
)::Int16
    material_id_new = material_id
    if material_id in felsic_cc_fertile ||
       material_id in felsic_cc_refrac ||
       material_id in solidified_granite ||
       material_id in solidified_rhyolite
        material_id_new = felsic_cc_pmolten[1]
    elseif material_id in mafic_cc_fertile ||
           material_id in mafic_cc_pmolten ||
           material_id in mafic_cc_refrac
        material_id_new = mafic_cc_pmolten[1]
    elseif material_id in ultramafic_fertile ||
           material_id in ultramafic_refrac
        material_id_new = ultramafic_pmolten[1]
    elseif material_id in solidified_gabbro ||
           material_id in solidified_basalt
        material_id_new = gabbro_pmolten[1]
    elseif material_id in solidified_layered_gabbro
        material_id_new = solidified_layered_gabbro_pmolten[1]
    end
    return material_id_new
end

""" Transform material id based on melt fraction.

This is used to visualize zones of partial melt.

# Arguments
- `material_id::Int16`: Material ID of the marker
- `felsic_cc_fertile::Vector{Int16}`: Array of felsic continental crust fertile IDs
- `felsic_cc_pmolten::Vector{Int16}`: Array of felsic continental crust partially 
  molten IDs
- `mafic_cc_fertile::Vector{Int16}`: Array of mafic continental crust fertile IDs
- `mafic_cc_pmolten::Vector{Int16}`: Array of mafic continental crust partially 
  molten IDs
- `ultramafic_fertile::Vector{Int16}`: Array of ultramafic mantle fertile IDs
- `ultramafic_pmolten::Vector{Int16}`: Array of ultramafic mantle partially 
  molten IDs
- `solidified_gabbro::Vector{Int16}`: Array of solidified gabbro IDs
- `gabbro_pmolten::Vector{Int16}`: Array of gabbro partially molten IDs
- `solidified_layered_gabbro::Vector{Int16}`: Array of solidified layered gabbro IDs
- `solidified_layered_gabbro_pmolten::Vector{Int16}`: Array of solidified layered 
  gabbro partially molten IDs

# Returns
- `Int16`: New material ID
"""
@inline function transform_partially_molten_marker_to_zero_melt_state(
    material_id::Int16,
    felsic_cc_fertile::Vector{Int16},
    felsic_cc_pmolten::Vector{Int16},
    mafic_cc_fertile::Vector{Int16},
    mafic_cc_pmolten::Vector{Int16},
    ultramafic_fertile::Vector{Int16},
    ultramafic_pmolten::Vector{Int16},
    solidified_gabbro::Vector{Int16},
    gabbro_pmolten::Vector{Int16},
    solidified_layered_gabbro::Vector{Int16},
    solidified_layered_gabbro_pmolten::Vector{Int16}
)::Int16
    material_id_new = material_id
    if material_id in felsic_cc_pmolten
        material_id_new = felsic_cc_fertile[1]
    elseif material_id in mafic_cc_pmolten
        material_id_new = mafic_cc_fertile[1]
    elseif material_id in ultramafic_pmolten
        material_id_new = ultramafic_fertile[1]
    elseif material_id in gabbro_pmolten
        material_id_new = solidified_gabbro[1]
    elseif material_id in solidified_layered_gabbro_pmolten
        material_id_new = solidified_layered_gabbro[1]
    end
    return material_id_new
end

end # module 