module MarkerTransformation

import EarthBox.ModelDataContainer: ModelData

@inline function transform_marker_to_air!(
    model::ModelData,
    imarker::Int,
    age_ma::Float64,
    matids_sticky_air::Vector{Int16}
)::Nothing
    matid_air = matids_sticky_air[1]
    basic_marker_reset!(model, imarker, age_ma, matid_air)
    return nothing
end

@inline function transform_marker_to_water!(
    model::ModelData,
    imarker::Int,
    age_ma::Float64,
    matids_sticky_water::Vector{Int16}
)::Nothing
    matid_water = matids_sticky_water[1]
    basic_marker_reset!(model, imarker, age_ma, matid_water)
    return nothing
end

@inline function transform_marker_to_sediment!(
    model::ModelData,
    imarker::Int,
    age_ma::Float64,
    matids_sediment::Vector{Int16}
)::Nothing
    matid_sediment = matids_sediment[1]
    basic_marker_reset!(model, imarker, age_ma, matid_sediment)
    return nothing
end

@inline function basic_marker_reset!(
    model::ModelData,
    imarker::Int,
    age_ma::Float64,
    new_material_id::Int16
)::Nothing
    @inbounds begin
        model.markers.arrays.material.marker_matid.array[imarker] = Int16(new_material_id)
        # Update friction coefficient
        model.markers.arrays.rheology.marker_fric.array[imarker] = 
            model.materials.arrays.mat_plastic.array[new_material_id, 3]
        # Update strain rate ratio
        model.markers.arrays.strain.marker_sr_ratio.array[imarker] = 1.0
        # Set strain to zero
        model.markers.arrays.strain.marker_GII.array[imarker] = 0.0
        model.markers.arrays.strain.marker_strain_plastic.array[imarker] = 0.0
        model.markers.arrays.strain.marker_strain_rate_plastic.array[imarker] = 0.0
        # Reset plastic failure flag
        model.markers.arrays.rheology.marker_pfailure.array[imarker] = 0.0f0
        # Reset maximum burial depth
        model.markers.arrays.material.marker_max_burial_depth.array[imarker] = 0.0
        # Reset serpentinization ratio
        model.markers.arrays.material.marker_serpentinization.array[imarker] = 0.0
        # Update age
        model.markers.arrays.strat.marker_age.array[imarker] = age_ma
    end
    return nothing
end

""" Transform marker to sediment or volcanics.

The distinction between sediment and volcanics is based on the thickness
of extruded material.

     ----> + x
     |
     v
    +y
            ------------ y_topo      -
               Volcanics             | Extrusive thickness
            ------------ y_sediment  -
               Sediments
"""
@inline function transform_marker_to_sediment_or_volcanics!(
    model::ModelData,
    imarker::Int,
    age_ma::Float64,
    y_topo::Float64,
    y_marker::Float64,
    iuse_extrusion::Int,
    extrusion_thickness::Float64,
    sedimentation_flag::Int,
    matids_sediment::Vector{Int16},
    matids_volcanics::Vector{Int16}
)::Nothing
    if iuse_extrusion == 0
        if sedimentation_flag == 1
            transform_marker_to_sediment!(model, imarker, age_ma, matids_sediment)
        end
    else
        y_sediment = y_topo + extrusion_thickness
        if y_marker > y_sediment
            if sedimentation_flag == 1
                transform_marker_to_sediment!(model, imarker, age_ma, matids_sediment)
            end
        else
            transform_marker_to_volcanics!(model, imarker, age_ma, matids_volcanics)
        end
    end
    return nothing
end

@inline function transform_marker_to_volcanics!(
    model::ModelData,
    imarker::Int,
    age_ma::Float64,
    matids_volcanics::Vector{Int16}
):Nothing
    @inbounds begin
        matid_lava = matids_volcanics[1]
        # Change marker type to volcanics
        model.markers.arrays.material.marker_matid.array[imarker] = matid_lava
        # Update friction coefficient
        model.markers.arrays.rheology.marker_fric.array[imarker] = 
            model.materials.arrays.mat_plastic.array[matid_lava, 3]
        # Reset strain rate ratio
        model.markers.arrays.strain.marker_sr_ratio.array[imarker] = 1.0
        # Reset strain
        model.markers.arrays.strain.marker_GII.array[imarker] = 0.0
        model.markers.arrays.strain.marker_strain_plastic.array[imarker] = 0.0
        model.markers.arrays.strain.marker_strain_rate_plastic.array[imarker] = 0.0
        # Reset normal stress
        model.markers.arrays.stress.marker_sxx.array[imarker] = 0.0
        # Reset shear stress
        model.markers.arrays.stress.marker_sxy.array[imarker] = 0.0
        # Reset viscosity
        model.markers.arrays.rheology.marker_eta.array[imarker] = 0.0
        # Reset normal strain rate
        model.markers.arrays.strain.marker_exx.array[imarker] = 0.0
        # Reset shear strain rate
        model.markers.arrays.strain.marker_exy.array[imarker] = 0.0
        # Update age
        model.markers.arrays.strat.marker_age.array[imarker] = age_ma
        # Reset plastic failure flag
        model.markers.arrays.rheology.marker_pfailure.array[imarker] = 0.0f0
        # Melt fraction; start at 1 since this is magma
        model.markers.arrays.melt.marker_meltfrac.array[imarker] = 1.0
        # This is the total extracted melt, set to zero
        model.markers.arrays.melt.marker_extracted_meltfrac.array[imarker] = 0.0
        # Reset incremental melt extraction
        model.markers.arrays.melt.marker_extractable_meltfrac.array[imarker] = 0.0
        # Reset serpentinization ratio
        model.markers.arrays.material.marker_serpentinization.array[imarker] = 0.0
        # Set temperature to emplacement temperature
        temperature_top = model.bcs.parameters.temperature.temperature_top.value
        model.markers.arrays.thermal.marker_TK.array[imarker] = temperature_top
    end
    return nothing
end

@inline function transform_rocks_to_water_or_air!(
    model::ModelData,
    y_sealevel::Float64,
    imarker::Int,
    age_ma::Float64,
    matids_sticky_air::Vector{Int16},
    matids_sticky_water::Vector{Int16}
)::Nothing
    @inbounds y_marker = model.markers.arrays.location.marker_y.array[imarker]
    if y_marker < y_sealevel
        transform_marker_to_air!(model, imarker, age_ma, matids_sticky_air)
    else
        transform_marker_to_water!(model, imarker, age_ma, matids_sticky_water)
    end
    return nothing
end

@inline function transform_marker_to_serpentinite!(
    model::ModelData,
    imarker::Int,
    age_ma::Float64,
    matid_serpentinite::Int
):Nothing
    @inbounds begin
        model.markers.arrays.material.marker_matid.array[imarker] = matid_serpentinite
        # Update friction coefficient
        model.markers.arrays.rheology.marker_fric.array[imarker] = 
            model.materials.arrays.mat_plastic.array[matid_serpentinite, 3]
        # Reset strain rate ratio
        model.markers.arrays.strain.marker_sr_ratio.array[imarker] = 1.0
        # Reset strain
        model.markers.arrays.strain.marker_GII.array[imarker] = 0.0
        model.markers.arrays.strain.marker_strain_plastic.array[imarker] = 0.0
        # Reset normal stress
        model.markers.arrays.stress.marker_sxx.array[imarker] = 0.0
        # Reset shear stress
        model.markers.arrays.stress.marker_sxy.array[imarker] = 0.0
        # Reset viscosity
        model.markers.arrays.rheology.marker_eta.array[imarker] = 0.0
        # Reset normal strain rate
        model.markers.arrays.strain.marker_exx.array[imarker] = 0.0
        # Reset shear strain rate
        model.markers.arrays.strain.marker_exy.array[imarker] = 0.0
        # Update age
        model.markers.arrays.strat.marker_age.array[imarker] = age_ma
        # Reset plastic failure flag
        model.markers.arrays.rheology.marker_pfailure.array[imarker] = 0.0f0
        # Melt fraction
        model.markers.arrays.melt.marker_meltfrac.array[imarker] = 0.0
        # This is the total extracted melt, set to zero
        model.markers.arrays.melt.marker_extracted_meltfrac.array[imarker] = 0.0
        # Reset incremental melt extraction
        model.markers.arrays.melt.marker_extractable_meltfrac.array[imarker] = 0.0
    end
    return nothing
end

@inline function transform_marker_to_salt!(
    model::ModelData,
    imarker::Int,
    age_ma::Float64,
    matids_salt::Vector{Int16}
):Nothing
    transform_marker_to_sediment!(model, imarker, age_ma, matids_salt)
    return nothing
end

end # module 