module MarkersStressRotation

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs

"""
    update!(model::ModelData, markers_spin::Vector{Float64})

Rotate marker stress using spin.

# Updated arrays from group `markers.arrays.stress`
- `marker_sxx::Vector{Float64}`: Marker normal stress (Pa)
- `marker_sxy::Vector{Float64}`: Marker shear stress (Pa)
"""
function update!(
    model::ModelData,
    inside_flags::Vector{Int8},
    markers_spin::Vector{Float64}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    marker_sxx = model.markers.arrays.stress.marker_sxx.array
    marker_sxy = model.markers.arrays.stress.marker_sxy.array
    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                spin = markers_spin[imarker]
                sxy = marker_sxy[imarker]
                sxx = marker_sxx[imarker]
                marker_sxy[imarker], marker_sxx[imarker] = rotate_marker_stress(
                    spin, timestep, sxy, sxx)
            end
        end
    end
    return nothing
end

""" Rotate marker stress using marker spin.

Rotation magnitude is calculated from spin rate:

    spin = 1/2(dvy/dx-dvx/dy)

Sign of rotation is positive for clockwise rotation when x axis is directed
rightward and y axis is directed downward.

New shear stress is calculated using the following equation:

sxy_new = 0.5(sxx_old - syy_old)*sin(2*spin*dt) + sxy_old*cos(2*spin*dt)

where sxx_old - syy_old = 2sxx_old.

New normal stress is calculated using the following equation:

sxx_new = sxx_old*(cos(spin*dt))^2 + syy_old*(sin(spin*dt))^2 - sxy_old*sin(2*spin*dt)

where sxx_old = -syy_old.
"""
@inline function rotate_marker_stress(
    spin::Float64,
    timestep::Float64,
    sxy_old::Float64,
    sxx_old::Float64
)::Tuple{Float64, Float64}
    spin_term1 = 2.0 * spin * timestep
    sxy_new = sxx_old * sin(spin_term1) + sxy_old * cos(spin_term1)
    spin_term2 = spin * timestep
    sxx_new = sxx_old * (cos(spin_term2)^2 - sin(spin_term2)^2) - sxy_old * sin(spin_term1)
    return sxy_new, sxx_new
end

end # module 