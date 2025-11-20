module GridTemperatureToMarkers

import EarthBox.ModelDataContainer: ModelData
import ..GridToMarker: get_marker_value

"""
    interpolate!(model::ModelData)::Nothing

Interpolate grid temperatures to markers.

# Arguments
- `model::ModelData`: Model data container

# Updated Arrays
## Updates from group `model.markers.arrays.thermal`:
- `marker_TK.array::Vector{Float64}` (marknum): Marker temperature in Kelvin
"""
function interpolate!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    # Get marker arrays
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array
    marknum = model.markers.parameters.distribution.marknum.value
    marker_temperature = model.markers.arrays.thermal.marker_TK.array

    # Get grid temperature array
    grid_temperature_tk1 = model.heat_equation.arrays.temperature.tk1.array

    # Interpolate temperatures to markers
    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                marker_temperature[imarker] = get_marker_value(
                    marker_yn[imarker],
                    marker_xn[imarker],
                    marker_dy[imarker],
                    marker_dx[imarker],
                    grid_temperature_tk1
                )
            end
        end
    end
    return nothing
end

end # module 