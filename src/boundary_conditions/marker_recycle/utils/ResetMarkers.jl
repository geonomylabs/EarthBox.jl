module ResetMarkers

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Markers.MarkerFriction.FrictionRandomizer: randomize_initial_friction_coefficient
import ..RandomMarkerArray: get_random_marker_array
import ..Sticky: get_sticky_material_ids
import ..Coordinates: MarkerRecycleArrays

struct SubsurfaceProps
    matids_subsurface_recycle::Vector{Int16}
    pressure_subsurface_recycle::Vector{Float64}
    temperature_subsurface_recycle::Vector{Float64}
end

function make_recycle_property_arrays(
    nrecycle_total::Int64
)::Tuple{Vector{Int16}, Vector{Float64}, Vector{Float64}}
    matids_recycle = fill(Int16(-1), nrecycle_total)
    pressure_recycle = zeros(Float64, nrecycle_total)
    temperature_recycle = zeros(Float64, nrecycle_total)
    return matids_recycle, pressure_recycle, temperature_recycle
end

function reset_recycled_marker_properties!(
    model::ModelData,
    marker_recycling::MarkerRecycleArrays,
    subsurface_props::SubsurfaceProps
)::Nothing
    reset_sticky_markers!(model, marker_recycling)
    reset_subsurface_markers!(model, marker_recycling, subsurface_props)
    return nothing
end

function reset_sticky_markers!(
    model::ModelData,
    marker_recycling::MarkerRecycleArrays
)::Nothing
    recycle_indices = marker_recycling.recycle_indices
    recycle_flags = marker_recycling.recycle_flags

    temperature_reset_sticky = model.bcs.parameters.temperature.temperature_top.value
    nrecycle_total = length(recycle_indices)
    matids_sticky = get_sticky_material_ids(model)
    matid_sticky_air = matids_sticky[1]
    
    for irecycle in 1:nrecycle_total
        imarker = recycle_indices[irecycle]
        recycle_flag = recycle_flags[irecycle]
        if recycle_flag == 0
            reset_marker!(
                model, imarker, matid_sticky_air,
                temperature_reset=temperature_reset_sticky
            )
        end
    end
    return nothing
end

function reset_subsurface_markers!(
    model::ModelData,
    marker_recycling::MarkerRecycleArrays,
    subsurface_props::SubsurfaceProps
)::Nothing
    matids_recycle = subsurface_props.matids_subsurface_recycle
    pressure_recycle = subsurface_props.pressure_subsurface_recycle
    temperature_recycle = subsurface_props.temperature_subsurface_recycle

    recycle_indices = marker_recycling.recycle_indices
    recycle_flags = marker_recycling.recycle_flags

    nrecycle_total = length(recycle_indices)
    marker_random = get_random_marker_array(model)
    iuse_random_fric, delta_fric_coef = get_friction_data(model)

    (
        _,
        stress_reset_subsurface,
        strain_rate_reset_subsurface,
        viscosity_reset_subsurface
    ) = get_default_reset_values(model)

    for irecycle in 1:nrecycle_total
        imarker = recycle_indices[irecycle]
        recycle_flag = recycle_flags[irecycle]
        if recycle_flag == 1
            matid_new = matids_recycle[irecycle]
            pressure_reset = pressure_recycle[irecycle]
            temperature_reset = temperature_recycle[irecycle]
            random_number = marker_random[imarker]
            reset_marker!(
                model, imarker, matid_new,
                temperature_reset=temperature_reset,
                pressure_reset=pressure_reset,
                strain_rate_reset=strain_rate_reset_subsurface,
                stress_reset=stress_reset_subsurface,
                viscosity_reset=viscosity_reset_subsurface,
                delta_fric_coef=delta_fric_coef,
                iuse_random_fric=iuse_random_fric,
                random_number=random_number
            )
        end
    end
    return nothing
end

function reset_marker!(
    model::ModelData,
    imarker::Int64,
    matid_reset::Int16;
    temperature_reset::Float64=0.0,
    pressure_reset::Float64=0.0,
    stress_reset::Float64=0.0,
    strain_rate_reset::Float64=0.0,
    viscosity_reset::Float64=0.0,
    delta_fric_coef::Float64=0.0,
    iuse_random_fric::Int64=0,
    random_number::Float64=0.0
)::Nothing
    reset_material_id!(model, imarker, matid_reset)
    reset_temperature!(model, imarker, temperature_reset)
    reset_pressure!(model, imarker, pressure_reset)
    reset_stress!(model, imarker, stress_reset)
    reset_strain_rate!(model, imarker, strain_rate_reset)
    reset_strain_rate_ratio!(model, imarker)
    reset_strain!(model, imarker)
    reset_melt!(model, imarker)
    reset_stratigraphy!(model, imarker)
    reset_rheology!(model, imarker, viscosity_reset)
    reset_friction_coefficients!(
        model, imarker, matid_reset, iuse_random_fric,
        delta_fric_coef, random_number
    )
    return nothing
end

function reset_material_id!(
    model::ModelData,
    imarker::Int64,
    matid_reset::Int16
)::Nothing
    marker_arrays = model.markers.arrays
    marker_arrays.material.marker_matid.array[imarker] = Int16(matid_reset)
    marker_arrays.material.marker_serpentinization.array[imarker] = 0.0
    return nothing
end

function reset_temperature!(
    model::ModelData,
    imarker::Int64,
    temperature_reset::Float64
)::Nothing
    marker_arrays = model.markers.arrays
    marker_arrays.thermal.marker_TK.array[imarker] = temperature_reset
    return nothing
end

function reset_pressure!(
    model::ModelData,
    imarker::Int64,
    pressure_reset::Float64
)::Nothing
    marker_arrays = model.markers.arrays
    marker_arrays.pressure.marker_pr.array[imarker] = pressure_reset
    return nothing
end

function reset_stress!(
    model::ModelData,
    imarker::Int64,
    stress_reset::Float64
)::Nothing
    marker_arrays = model.markers.arrays
    stress = marker_arrays.stress
    stress.marker_sxx.array[imarker] = stress_reset
    stress.marker_sxy.array[imarker] = stress_reset
    return nothing
end

function reset_strain_rate!(
    model::ModelData,
    imarker::Int64,
    strain_rate_reset::Float64
)::Nothing
    marker_arrays = model.markers.arrays
    strain = marker_arrays.strain
    strain.marker_exx.array[imarker] = strain_rate_reset
    strain.marker_exy.array[imarker] = strain_rate_reset
    return nothing
end

function reset_strain_rate_ratio!(
    model::ModelData,
    imarker::Int64
)::Nothing
    marker_arrays = model.markers.arrays
    strain = marker_arrays.strain
    strain.marker_sr_ratio.array[imarker] = 0.0
    return nothing
end

function reset_strain!(
    model::ModelData,
    imarker::Int64
)::Nothing
    marker_arrays = model.markers.arrays
    strain = marker_arrays.strain
    strain.marker_GII.array[imarker] = 0.0
    strain.marker_strain_plastic.array[imarker] = 0.0
    strain.marker_strain_rate_plastic.array[imarker] = 0.0
    return nothing
end

function reset_melt!(
    model::ModelData,
    imarker::Int64
)::Nothing
    marker_arrays = model.markers.arrays
    melt = marker_arrays.melt
    melt.marker_meltfrac.array[imarker] = 0.0
    melt.marker_extracted_meltfrac.array[imarker] = 0.0
    melt.marker_extractable_meltfrac.array[imarker] = 0.0
    return nothing
end

function reset_stratigraphy!(
    model::ModelData,
    imarker::Int64
)::Nothing
    marker_arrays = model.markers.arrays
    marker_arrays.strat.marker_age.array[imarker] = 0.0
    return nothing
end

function reset_rheology!(
    model::ModelData,
    imarker::Int64,
    viscosity_reset::Float64
)::Nothing
    marker_arrays = model.markers.arrays
    rheology = marker_arrays.rheology
    rheology.marker_eta.array[imarker] = viscosity_reset
    rheology.marker_pfailure.array[imarker] = 0.0f0
    return nothing
end

function reset_friction_coefficients!(
    model::ModelData,
    imarker::Int64,
    matid_reset::Int16,
    iuse_random_fric::Int64,
    delta_fric_coef::Float64,
    random_number::Float64
)::Nothing
    marker_arrays = model.markers.arrays
    rheology = marker_arrays.rheology

    mat_plastic = model.materials.arrays.mat_plastic.array
    friction_coefficient = mat_plastic[matid_reset,3]

    if iuse_random_fric == 1
        friction_coefficient = randomize_initial_friction_coefficient(
            friction_coefficient, delta_fric_coef, random_number)
    end

    rheology.marker_fric_ini.array[imarker] = friction_coefficient
    rheology.marker_fric.array[imarker] = friction_coefficient
    return nothing
end

function get_friction_data(model::ModelData)::Tuple{Int, Float64}
    random_friction = model.materials.parameters.random_friction
    iuse_random_fric = random_friction.iuse_random_fric.value
    delta_fric_coef = random_friction.delta_fric_coef.value
    return (iuse_random_fric, delta_fric_coef)
end

function get_default_reset_values(
    model::ModelData
)::Tuple{Float64, Float64, Float64, Float64}
    temperature_reset_sticky = model.bcs.parameters.temperature.temperature_top.value
    stress_reset = 1.96e6
    strain_rate_reset = 1e-14
    viscosity_reset = 1e20
    
    return (
        temperature_reset_sticky,
        stress_reset,
        strain_rate_reset,
        viscosity_reset
    )
end

end # module ResetMarkers 