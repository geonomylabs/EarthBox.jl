module PlumeUpdate

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConversionFuncs: seconds_to_years
import ....MarkerMaterials.InitManager.Plume: in_plume_region_check
import ....MarkerMaterials.InitManager.InitStructs: PlumeGeometry

""" Inject plume if plume start time is reached.
"""
function inject_plume!(model::ModelData)
    _inject_plume = check_for_plume_injection(model)
    if _inject_plume
        update_temperature_for_plume!(model)
    end
end

""" Check for plume injection and update plume state.
"""
function check_for_plume_injection(model::ModelData)
    plume_start_time = model.geometry.parameters.plume.plume_start_time.value
    iplume_state = model.geometry.parameters.plume.iplume_state.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_myr = seconds_to_years(timesum)/1e6
    _inject_plume = false
    if iplume_state == 0 && timesum_myr >= plume_start_time
        _inject_plume = true
        model.geometry.parameters.plume.iplume_state.value = 1
    end
    return _inject_plume
end

""" Update temperature for plume.

Updated arrays:
- `model.markers.arrays.thermal.marker_TK.array`: Temperature of markers in Kelvin.
"""
function update_temperature_for_plume!(model::ModelData)
    marknum                 = model.markers.parameters.distribution.marknum.value
    iuse_plume              = model.geometry.parameters.plume.iuse_plume.value
    delta_temperature_plume = model.geometry.parameters.plume.delta_temperature_plume.value
    plume_radius            = model.geometry.parameters.plume.plume_radius.value
    plume_center_x          = model.geometry.parameters.plume.plume_center_x.value
    plume_center_y          = model.geometry.parameters.plume.plume_center_y.value
    plume_head_thick        = model.geometry.parameters.plume.plume_head_thick.value
    
    marker_x  = model.markers.arrays.location.marker_x.array
    marker_y  = model.markers.arrays.location.marker_y.array
    marker_TK = model.markers.arrays.thermal.marker_TK.array

    if iuse_plume == 1
        plume_geometry = PlumeGeometry(
            plume_radius = plume_radius,
            plume_center_x = plume_center_x,
            plume_center_y = plume_center_y,
            plume_head_thick = plume_head_thick
        )        
        Threads.@threads for imarker in 1:marknum
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            in_plume = in_plume_region_check(x_marker, y_marker, plume_geometry)
            if in_plume
                marker_TK[imarker] = marker_TK[imarker] + delta_temperature_plume
            end
        end
    end
end

end # module 