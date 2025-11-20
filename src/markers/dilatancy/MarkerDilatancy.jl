module MarkerDilatancy

using Printf
using EarthBox.ModelDataContainer: ModelData

"""
    initialize!(model::ModelData)::Nothing

Initialize marker dilatancy for plasticity models. This initializes the `marker_dilatation_angle`
array in the model data container and is required for all models that use dilatancy.

# Arguments
- `model::ModelData`: The model data container containing the model parameters and arrays.

# Returns
- `Nothing`
"""
function initialize!(model::ModelData)::Nothing
    initialize_marker_dilatancy!(model)
    nothing
end

function initialize_marker_dilatancy!(model::ModelData)::Nothing
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_dilatation_angle = model.markers.arrays.rheology.marker_dilatation_angle.array
    mat_plastic = model.materials.arrays.mat_plastic.array
    marknum = model.markers.parameters.distribution.marknum.value

    for imarker in 1:marknum
        matid = marker_matid[imarker]
        marker_dilatation_angle[imarker] = Float32(mat_plastic[matid,7])
    end

    return nothing
end

end # module 