module FractureZoneProperties

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Markers.MarkerMaterials.InitManager.FractureZone: FractureZoneGeometry
import EarthBox.Markers.MarkerMaterials.InitManager.FractureZone: FractureZoneMaterialIDs
import EarthBox.Markers.MarkerMaterials.InitManager.FractureZone: define_geometry as define_fz_geometry
import EarthBox.Markers.MarkerMaterials.InitManager.FractureZone: define_material_ids as define_fz_material_ids
import EarthBox.Markers.MarkerMaterials.InitManager.FractureZone: define_material as define_fz_material
import EarthBox.Markers.MarkerTemperature.InitManager.FractureZone: ThermalParameters
import EarthBox.Markers.MarkerTemperature.InitManager.FractureZone: get_thermal_parameters as get_fz_thermal_parameters
import EarthBox.Markers.MarkerTemperature.InitManager.FractureZone: calculate_temperature as calculate_fz_temperature
import ..Coordinates: MarkerRecycleArrays
import ..PressureReset: calculate_reset_pressure
import ..ResetMarkers: SubsurfaceProps

function get_subsurface_reset_properties_fracture_zone(
    model::ModelData,
    marker_recycling::MarkerRecycleArrays
)::SubsurfaceProps
    recycle_indices = marker_recycling.recycle_indices
    recycle_flags = marker_recycling.recycle_flags
    nrecycle_total = length(recycle_indices)

    matids_recycle = fill(-1, nrecycle_total)
    pressure_recycle = zeros(nrecycle_total)
    temperature_recycle = zeros(nrecycle_total)

    marker_arrays = model.markers.arrays
    location = marker_arrays.location
    marker_x = location.marker_x.array
    marker_y = location.marker_y.array

    (
        fz_material_ids, fz_geometry, thermal_parameters
    ) = get_fracture_zone_model_info(model)

    Threads.@threads for irecycle in 1:nrecycle_total
        recycle_flag = recycle_flags[irecycle]
        if recycle_flag == 1
            imarker = recycle_indices[irecycle]
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            matids_recycle[irecycle] = define_fz_material(
                x_marker, y_marker, fz_material_ids, fz_geometry)
            pressure_recycle[irecycle] = calculate_reset_pressure(y_marker)
            temperature_recycle[irecycle] = calculate_fz_temperature(
                y_marker, x_marker, thermal_parameters)
        end
    end

    return SubsurfaceProps(
        matids_recycle,
        pressure_recycle,
        temperature_recycle
    )
end

function get_fracture_zone_model_info(model::ModelData
)::Tuple{FractureZoneMaterialIDs, FractureZoneGeometry, ThermalParameters}
    fz_material_ids = define_fz_material_ids(model)
    fz_geometry = define_fz_geometry(model)
    thermal_parameters = get_fz_thermal_parameters(model)
    return (fz_material_ids, fz_geometry, thermal_parameters)
end

end # module 