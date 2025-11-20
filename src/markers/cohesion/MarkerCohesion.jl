module MarkerCohesion

using EarthBox.ModelDataContainer: ModelData

"""
    initialize!(model::ModelData)::Nothing

Initialize marker cohesion for plasticity models. This initializes the `marker_cohesion`
array in the model data container and is required for all models that use plasticity.

# Arguments
- `model::ModelData`: The model data container containing the model parameters and arrays.

# Returns
- `Nothing`
"""
function initialize!(model::ModelData)::Nothing
    initialize_marker_cohesion!(model)
    nothing
end

"""
    initialize_marker_cohesion!(model::ModelData)

Initialize marker cohesion.

# Arguments
- `model::ModelData`: The model containing marker and material data

# Updates
- `model.markers.arrays.rheology.marker_cohesion.array`: Array of marker cohesion values (Pa)
"""
function initialize_marker_cohesion!(model::ModelData)
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_cohesion = model.markers.arrays.rheology.marker_cohesion.array
    mat_plastic = model.materials.arrays.mat_plastic.array
    marknum = model.markers.parameters.distribution.marknum.value
    for imarker in 1:marknum
        matid = marker_matid[imarker]
        marker_cohesion[imarker] = mat_plastic[matid][1]
    end
end

end # module 