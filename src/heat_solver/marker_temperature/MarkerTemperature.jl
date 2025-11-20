module MarkerTemperature

import EarthBox.ModelDataContainer: ModelData
import EarthBox.GridFuncs: check_in_domain
import EarthBox.Interpolation: GridToMarker
import EarthBox.BoundaryConditions.MarkerRecycle.Sticky: get_sticky_material_ids

""" Add remaining grid temperature changes to marker temperature.

The remaining conductive grid temperature changes on the marker 
delta_grid_temperature_remaining is obtained by interpolating the remaining
conductive grid temperature changes `dtk1` on the basic grid to the marker
location using bilinear interpolation.

# Updated Arrays
- `marker_TK`: Array of marker temperatures (K)
"""
function add_remaining_grid_change_to_marker_temperature!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    marknum            = model.markers.parameters.distribution.marknum.value
    marker_temperature = model.markers.arrays.thermal.marker_TK.array

    dtk1      = model.heat_equation.arrays.temperature.dtk1.array
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                dx_upr_left = marker_dx[imarker]
                dy_upr_left = marker_dy[imarker]
                ix_upr_left = marker_xn[imarker]
                iy_upr_left = marker_yn[imarker]
            end
            delta_grid_temperature_remaining = 
                GridToMarker.get_marker_value_from_basic_grid_array(
                    iy_upr_left, ix_upr_left, dy_upr_left, dx_upr_left, dtk1)
            @inbounds marker_temperature[imarker] += delta_grid_temperature_remaining
        end
    end
    return nothing
end

""" Interpolate temperature change on basic grid to markers.

# Updated Arrays
- `marker_TK`: Array of marker temperatures (K)
"""
function apply_sticky_temperature_correction!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    temperature_top = model.bcs.parameters.temperature.temperature_top.value
    marker_temperature = model.markers.arrays.thermal.marker_TK.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    matids_sticky = get_sticky_material_ids(model)

    for imarker in 1:marknum
        if inside_flags[imarker] == 1
            temperature_kelvin = marker_temperature[imarker]
            matid = marker_matid[imarker]
            if matid in matids_sticky && temperature_kelvin != temperature_top
                marker_temperature[imarker] = temperature_top
            end
        end
    end
end

end # module 