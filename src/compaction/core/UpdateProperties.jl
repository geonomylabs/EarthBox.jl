module UpdateProperties

import EarthBox.ModelDataContainer: ModelData

""" Initialize marker compaction properties in model data.

This function transfers compaction properties in material arrays to marker arrays.

# Updated Arrays
- `model.markers.arrays.material.marker_porosity_initial`: Initial porosity of markers at the 
  sediment-water interface used in the exponential decay model.
- `model.markers.arrays.material.marker_decay_depth`: Porosity decay depth of markers used in 
  the exponential decay model.
"""
function reset_compaction_properties!(model::ModelData)
    marker_porosity_initial = model.markers.arrays.material.marker_porosity_initial.array
    marker_decay_depth = model.markers.arrays.material.marker_decay_depth.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    mat_compaction = model.materials.arrays.mat_compaction.array
    nmarkers = size(marker_matid, 1)

    Threads.@threads for imarker in 1:nmarkers
        matid = marker_matid[imarker]
        porosity_initial = mat_compaction[matid, 1]
        decay_depth = mat_compaction[matid, 2]
        marker_porosity_initial[imarker] = porosity_initial
        marker_decay_depth[imarker] = decay_depth
    end
end

end