module RandomFrictionAngles

import EarthBox.ModelDataContainer: ModelData
using Printf
import Random

"""
    randomize_initial_friction_angles!(model::ModelData)::Nothing

Randomize initial marker friction angles by calling the loop function.

This randomizes the initial friction angle θ°ₘ for each marker that is stored as the sine of this angle
in the array `marker_fric_ini`.
"""
function randomize_initial_friction_angles!(model::ModelData)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    marker_random = Random.rand(marknum)
    randomize_initial_friction_angles_loop!(model, marker_random)
    return nothing
end

"""
    randomize_initial_friction_angles_loop!(
        model::ModelData,
        marker_random::Vector{Float64}
    )::Nothing

Randomize initial marker friction angles.

This randomizes the initial friction angle θ°ₘ for each marker and stores the sine of this angle
(i.e. friction coefficient) in the array `marker_fric_ini`. This function assumes that sticky air
material id's are 0 and 1, which are excluded from the randomization process.

# Background
The friction angle θₘ used in the yield stress equation is described with the following equation:

θₘ = θ°ₘ + (θᶠₘ - θ°ₘ)/(εᶠₘ - ε°ₘ)*εₘ

where m is the index of the marker, θ°ₘ and θᶠₘ are the initial and final friction angles,
respectively, εₘ is the current plastic strain and ε°ₘ and εᶠₘ are the initial and final
plastic strains, respectively.

The randomization of the initial friction angle is performed by adding a random perturbation
to the initial friction angle θ°ₘ for each marker as follows:

θ°ₘ = θ°ₘ + (0.5 - rₘ)*f

where rₘ is a random number between 0 and 1 for each marker and f is the user defined
randomization factor. A typical value of f is 10.
"""
function randomize_initial_friction_angles_loop!(
    model::ModelData,
    marker_random::Vector{Float64}
)::Nothing
    f_random = model.materials.parameters.random_friction.randomization_factor.value
    marknum = model.markers.parameters.distribution.marknum.value

    marker_matid = model.markers.arrays.material.marker_matid.array
    mat_plastic = model.materials.arrays.mat_plastic.array
    marker_fric_ini = model.markers.arrays.rheology.marker_fric_ini.array

    pi_term = π / 180.0
    inv_pi_term = 180.0 / π
    Threads.@threads for imarker in 1:marknum
        @inbounds begin
            mitype = marker_matid[imarker]
            rnum = marker_random[imarker]
            friction_coefficient_initial = mat_plastic[mitype][2]
        end
        friction_angle_initial = asin(friction_coefficient_initial) * inv_pi_term
        if 1 < mitype < 20
            friction_angle_initial_randomized = friction_angle_initial + (0.5 - rnum) * f_random
            @inbounds marker_fric_ini[imarker] = sin(friction_angle_initial_randomized * pi_term)
        end
    end
    
    return nothing
end

"""
    get_min_and_max_friction_angles(model::ModelData)::Tuple{Float64, Float64}

Calculate minimum and maximum friction angles excluding air and strong zones.
"""
function get_min_and_max_friction_angles(model::ModelData)::Tuple{Float64, Float64}
    marknum = model.markers.parameters.distribution.marknum.value
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_fric_ini = model.markers.arrays.rheology.marker_fric_ini.array
    
    theta_min = 1e34
    theta_max = -1e34
    
    for imarker in 1:marknum
        mitype = marker_matid[imarker]
        mfric = marker_fric_ini[imarker]
        theta = asin(mfric) * 180.0 / π
        
        # Avoid sticky air and water and strong zone
        if 1 < mitype < 20
            if theta < theta_min
                theta_min = theta
            end
            if theta > theta_max
                theta_max = theta
            end
        end
    end
    
    return theta_min, theta_max
end

end # module 