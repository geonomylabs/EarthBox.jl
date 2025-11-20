module MarkerPreexponential

import EarthBox.ModelDataContainer: ModelData

"""
    initialize!(model::ModelData)::Nothing

Initialize marker pre-exponential term for dislocation creep models. This 
initializes the `marker_preexp` array in the model data container and is required
for all models that use dislocation creep. The initialization of marker pre-exponential
terms allows strain softening to be applied to the pre-exponential term.

# Arguments
- `model::ModelData`: The model data container containing the model parameters and arrays.

# Returns
- `Nothing`

"""
function initialize!(model::ModelData)::Nothing
    initialize_marker_preexponential!(model)
    nothing
end

"""
    initialize_marker_preexponential!(model::ModelData)::Nothing

Initialize marker pre-exponential for dislocation creep.

# Updated Arrays
- `model.markers.arrays.rheology.marker_preexp.array`: 
    - Pre-exponential factor for dislocation creep of marker (1/s/MPa^n).
"""
function initialize_marker_preexponential!(model::ModelData)::Nothing
    marker_matid = model.markers.arrays.material.marker_matid.array
    mat_flow = model.materials.arrays.mat_flow.array
    marker_preexp = model.markers.arrays.rheology.marker_preexp.array
    marknum = model.markers.parameters.distribution.marknum.value
    
    for imarker in 1:marknum
        matid = marker_matid[imarker]
        marker_preexp[imarker] = mat_flow[matid, 2]
    end
    
    nothing
end

end # module