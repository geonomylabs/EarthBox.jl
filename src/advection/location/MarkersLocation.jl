module MarkersLocation
import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs

""" Displace markers based on their velocities.

# Updated arrays from group `markers.arrays.location`
- `marker_x::Vector{Float64}`: x-location of markers (m)
- `marker_y::Vector{Float64}`: y-location of markers (m)
"""
function update!(
    model::ModelData,
    inside_flags::Vector{Int8},
    markers_vx::Vector{Float64},
    markers_vy::Vector{Float64}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                x_marker = marker_x[imarker]
                y_marker = marker_y[imarker]
                vx = markers_vx[imarker]
                vy = markers_vy[imarker]
                marker_x[imarker] = advect_marker(x_marker, vx, timestep)
                marker_y[imarker] = advect_marker(y_marker, vy, timestep)
            end
        end
    end
    return nothing
end

@inline function advect_marker(
    marker_loc::Float64,
    velocity::Float64,
    timestep::Float64
)::Float64
    return marker_loc + timestep * velocity
end 

end # module