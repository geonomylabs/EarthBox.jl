module Regular

import EarthBox.ModelDataContainer: ModelData

"""
    initialize!(model::ModelData)

Update marker friction coefficients.

# Updated Arrays
## Updated arrays from group `model.markers.arrays.rheology`
- `marker_friction.array::Vector{Float64}`:
    - Array of friction coefficients (sine of friction angle) of markers.
- `marker_fric_ini.array::Vector{Float64}`:
    - Array of initial friction coefficients (sine of friction angle) of markers.
"""
function initialize!(model::ModelData)
    marker_matid = model.markers.arrays.material.marker_matid.array
    mat_plastic = model.materials.arrays.mat_plastic.array

    marker_fric_ini = model.markers.arrays.rheology.marker_fric_ini.array
    marker_fric = model.markers.arrays.rheology.marker_fric.array

    matids_sticky_air = model.materials.dicts.matid_types["StickyAir"]
    matids_sticky_water = model.materials.dicts.matid_types["StickyWater"]

    marknum = model.markers.parameters.distribution.marknum.value
    for imarker in 1:marknum
        matid = marker_matid[imarker]
        friction_coefficient = mat_plastic[matid, 3]  # Julia is 1-indexed
        # Avoid sticky air and water (This is unnecessary and can be removed)
        # This is a relict of the randomized friction option
        if matid in matids_sticky_air || matid in matids_sticky_water
            marker_fric_ini[imarker] = mat_plastic[matid, 3]  # Julia is 1-indexed
        else
            marker_fric_ini[imarker] = friction_coefficient
        end
        marker_fric[imarker] = marker_fric_ini[imarker]
    end
    
    return nothing
end

end # module 