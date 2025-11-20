module FailureProperties

using EarthBox.ModelDataContainer: ModelData
using EarthBox.GridFuncs: check_in_domain
using ..BoundaryFriction: apply_boundary_friction

function calculate_failure_properties!(
    model::ModelData,
    inside_flags::Vector{Int8},
    use_boundary_friction_plasticity_model_sandbox::Bool
)::Nothing
    calculate_failure_properties_loop!(model, inside_flags)
    if use_boundary_friction_plasticity_model_sandbox
        apply_boundary_friction(model)
    end
    return nothing
end

function calculate_failure_properties_loop!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    # Unpack data structures for better performance
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_strain_plastic = model.markers.arrays.strain.marker_strain_plastic.array

    marker_melt_damage = model.markers.arrays.strain.marker_melt_damage.array
    iuse_melt_damage = model.materials.parameters.melt_damage.iuse_melt_damage.value

    marker_fric_ini = model.markers.arrays.rheology.marker_fric_ini.array
    marker_cohesion = model.markers.arrays.rheology.marker_cohesion.array
    marker_fric = model.markers.arrays.rheology.marker_fric.array

    mat_plastic = model.materials.arrays.mat_plastic.array

    # Print min/max of melt damage factor if melt damage is enabled
    if iuse_melt_damage == 1
        min_melt_damage = minimum(marker_melt_damage)
        max_melt_damage = maximum(marker_melt_damage)
    end
    
    marknum = model.markers.parameters.distribution.marknum.value
    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            mat_id = marker_matid[imarker]
            strain = marker_strain_plastic[imarker]
            sine_friction_angle_initial_marker = marker_fric_ini[imarker]

            cohesion_initial = mat_plastic[mat_id, 1]
            cohesion_final = mat_plastic[mat_id, 2]
            sine_friction_angle_initial_mat = mat_plastic[mat_id, 3]
            sine_friction_angle_final_mat = mat_plastic[mat_id, 4]
            strain_initial = mat_plastic[mat_id, 5]
            strain_final = mat_plastic[mat_id, 6]

            sine_friction_angle_final = update_final_friction_angle(
                sine_friction_angle_initial_mat,
                sine_friction_angle_final_mat,
                sine_friction_angle_initial_marker
            )

            (cohesion, sine_friction_angle) = strain_weakening_or_hardening(
                cohesion_initial, cohesion_final,
                sine_friction_angle_initial_marker,
                sine_friction_angle_final,
                strain_initial, strain_final, strain
            )

            if iuse_melt_damage == 1
                melt_damage_factor = marker_melt_damage[imarker]
                sine_friction_angle = sine_friction_angle / melt_damage_factor
                cohesion = cohesion / melt_damage_factor
            end

            marker_cohesion[imarker] = cohesion
            marker_fric[imarker] = sine_friction_angle
        end
    end
    return nothing
end
 
function update_final_friction_angle(
    sine_friction_angle_initial_mat::Float64,
    sine_friction_angle_final_mat::Float64,
    sine_friction_angle_initial_marker::Float64
)::Float64
    if sine_friction_angle_initial_mat > 0.0 && sine_friction_angle_final_mat > 0.0
        fac = sine_friction_angle_final_mat / sine_friction_angle_initial_mat
        sine_friction_angle_final = sine_friction_angle_initial_marker * fac
    else
        sine_friction_angle_final = sine_friction_angle_final_mat
    end
    return sine_friction_angle_final
end

function strain_weakening_or_hardening(
    cohesion_initial::Float64,
    cohesion_final::Float64,
    sine_friction_angle_initial::Float64,
    sine_friction_angle_final::Float64,
    strain_initial::Float64,
    strain_final::Float64,
    strain::Float64
)::Tuple{Float64, Float64}
    # Checking yielding criterion for strain weakening/hardening
    cohesion = cohesion_initial
    sine_friction_angle = sine_friction_angle_initial
    if strain >= strain_final
        cohesion = cohesion_final
        sine_friction_angle = sine_friction_angle_final
    elseif strain_initial <= strain < strain_final
        cohesion = calc_cohesion(
            cohesion_initial,
            cohesion_final,
            strain_initial,
            strain_final,
            strain
        )
        sine_friction_angle = calc_friction(
            sine_friction_angle_initial,
            sine_friction_angle_final,
            strain_initial,
            strain_final,
            strain
        )
    end
    return (cohesion, sine_friction_angle)
end

function calc_cohesion(
    cohesion_initial::Float64,
    cohesion_final::Float64,
    strain_initial::Float64,
    strain_final::Float64,
    strain::Float64
)::Float64
    cohesion = (
        cohesion_initial
        + (cohesion_final - cohesion_initial) /
        (strain_final - strain_initial) * (strain - strain_initial)
    )
    return cohesion
end

function calc_friction(
    sine_friction_angle_initial::Float64,
    sine_friction_angle_final::Float64,
    strain_initial::Float64,
    strain_final::Float64,
    strain::Float64
)::Float64
    sine_friction_angle = (
        sine_friction_angle_initial
        + (sine_friction_angle_final - sine_friction_angle_initial) /
        (strain_final - strain_initial) * (strain - strain_initial)
    )
    return sine_friction_angle
end

end # module 