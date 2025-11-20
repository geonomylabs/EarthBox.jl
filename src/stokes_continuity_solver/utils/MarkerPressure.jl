module MarkerPressure

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Interpolation: MarkerGridMapping
import EarthBox.Interpolation: GridToMarker
import EarthBox: GridFuncs

""" Interpolate grid pressure to markers.

Updated Arrays
==============
model.markers.arrays.pressure
----------------------------
marker_pr: Array{Float64,1}
    - Marker pressure (Pa).
"""
function calculate_marker_pressure!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    xnum    = model.grids.parameters.geometry.xnum.value
    ynum    = model.grids.parameters.geometry.ynum.value
    marknum = model.markers.parameters.distribution.marknum.value

    marker_pr = model.markers.arrays.pressure.marker_pr.array
    marker_x  = model.markers.arrays.location.marker_x.array
    marker_y  = model.markers.arrays.location.marker_y.array
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    gridy_vx  = model.grids.arrays.staggered_vx.gridy_vx.array
    ystp_vx   = model.grids.arrays.staggered_vx.ystp_vx.array
    gridx_vy  = model.grids.arrays.staggered_vy.gridx_vy.array
    xstp_vy   = model.grids.arrays.staggered_vy.xstp_vy.array
    pr1       = model.stokes_continuity.arrays.pressure.pr1.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                x_marker = marker_x[imarker]
                y_marker = marker_y[imarker]
                yn = marker_yn[imarker]
                xn = marker_xn[imarker]
            end
            # Define indexes and normalized distances for upper left pressure
            # grid node in the cell that contains the marker. Note that
            # indices and distances for the basic grid have been pre-computed.
            yn_ul_pr = MarkerGridMapping.upr_left_index_pressure(ynum, yn, y_marker, gridy_vx)
            xn_ul_pr = MarkerGridMapping.upr_left_index_pressure(xnum, xn, x_marker, gridx_vy)
            dy_ul_pr = MarkerGridMapping.upr_left_dist_pressure(yn_ul_pr, y_marker, gridy_vx, ystp_vx)
            dx_ul_pr = MarkerGridMapping.upr_left_dist_pressure(xn_ul_pr, x_marker, gridx_vy, xstp_vy)
            # Interpolate grid pressure to marker
            @inbounds marker_pr[imarker] = GridToMarker.get_marker_value(
                yn_ul_pr, xn_ul_pr, dy_ul_pr, dx_ul_pr, pr1)
        end
    end
    return nothing
end

end # module 