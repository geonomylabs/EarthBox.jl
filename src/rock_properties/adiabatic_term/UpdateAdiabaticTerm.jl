module UpdateAdiabaticTerm

using EarthBox.ModelDataContainer.ModelDataContainer: ModelData

"""
    update_adiabatic_term!(model::ModelData)

Calculate rock properties for all markers.

# Arguments
- `model::ModelData`:
    - The model data container containing marker and material information

# Updated arrays in model.markers.arrays.thermal:
- `marker_ha.array`:
    - Adiabatic term (expansivity × temperature) for each marker.

"""
function update_adiabatic_term!(model::ModelData, inside_flags::Vector{Int8})
    marker_ha = model.markers.arrays.thermal.marker_ha.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                material_id = marker_matid[imarker]
                temperature_kelvins = model.markers.arrays.thermal.marker_TK.array[imarker]
                expansivity = model.materials.arrays.mat_rho.array[material_id, 2]  # Note: Julia uses 1-based indexing
            end
            adiabatic_term = expansivity * temperature_kelvins
            @inbounds marker_ha[imarker] = adiabatic_term
        end
    end
end

end # module 