module Randomized

using Random
import EarthBox.ModelDataContainer: ModelData
import ...FrictionRandomizer: randomize_initial_friction_coefficient

"""
    initialize!(model)

Update marker friction coefficients via randomization.

# Updated Arrays
## Updated arrays from group `model.markers.arrays.rheology`
- `marker_friction.array::Vector{Float64}`:
    - Array of friction coefficients (sine of friction angle) of markers.
- `marker_fric_ini.array::Vector{Float64}`:
    - Array of initial friction coefficients (sine of friction angle) of markers.
"""
function initialize!(model::ModelData)
    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    marker_random = rand(mxnum*mynum)

    marker_matid = model.markers.arrays.material.marker_matid.array
    mat_plastic = model.materials.arrays.mat_plastic.array
    delta_fric_coef = model.materials.parameters.random_friction.delta_fric_coef.value

    marker_fric_ini = model.markers.arrays.rheology.marker_fric_ini.array
    marker_fric = model.markers.arrays.rheology.marker_fric.array

    matids_sticky_air = model.materials.dicts.matid_types["StickyAir"]
    matids_sticky_water = model.materials.dicts.matid_types["StickyWater"]

    marknum = model.markers.parameters.distribution.marknum.value
    for imarker in 1:marknum
        matid = marker_matid[imarker]
        friction_coefficient = mat_plastic[matid][3]  # Julia is 1-indexed
        random_number = marker_random[imarker]
        friction_coefficient = randomize_initial_friction_coefficient(
            friction_coefficient, delta_fric_coef, random_number)
        
        # Avoid sticky air and water
        if matid in matids_sticky_air || matid in matids_sticky_water
            marker_fric_ini[imarker] = mat_plastic[matid][3]  # Julia is 1-indexed
        else
            marker_fric_ini[imarker] = friction_coefficient
        end
        marker_fric[imarker] = marker_fric_ini[imarker]
    end
    
    return nothing
end

end # module 