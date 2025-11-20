module MarkerStrainRate

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Interpolation: MarkerGridMapping
import EarthBox.Interpolation: GridToMarker
import EarthBox: GridFuncs

""" Interpolate normal and shear grid strain rates (exx, exy) to markers.

Normal and shear grid strain rates (exx and exy) were calculated from new
Stokes velocity solution.

Updated Arrays
==============
model.markers.arrays.strain
--------------------------
marker_exx: Array{Float64,1}
    Marker normal strain rate (1/s).

marker_exy: Array{Float64,1}
    Marker shear strain rate (1/s).
"""
function calculate_marker_strain_rate!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    xnum    = model.grids.parameters.geometry.xnum.value
    ynum    = model.grids.parameters.geometry.ynum.value
    marknum = model.markers.parameters.distribution.marknum.value

    marker_x   = model.markers.arrays.location.marker_x.array
    marker_y   = model.markers.arrays.location.marker_y.array
    marker_xn  = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn  = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_dx  = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy  = model.markers.arrays.grid_marker_relationship.marker_dy.array
    marker_exx = model.markers.arrays.strain.marker_exx.array
    marker_exy = model.markers.arrays.strain.marker_exy.array
    gridy_vx   = model.grids.arrays.staggered_vx.gridy_vx.array
    ystp_vx    = model.grids.arrays.staggered_vx.ystp_vx.array
    gridx_vy   = model.grids.arrays.staggered_vy.gridx_vy.array
    xstp_vy    = model.grids.arrays.staggered_vy.xstp_vy.array
    exx        = model.stokes_continuity.arrays.strain_rate_and_spin.exx.array
    exy        = model.stokes_continuity.arrays.strain_rate_and_spin.exy.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                x_marker = marker_x[imarker]
                y_marker = marker_y[imarker]
                yn = marker_yn[imarker]
                xn = marker_xn[imarker]
                dy = marker_dy[imarker]
                dx = marker_dx[imarker]
            end
            # Define indexes and normalized distances for upper left pressure
            # grid node in the cell that contains the marker. Notes that
            # indices and distances for the basic grid have been pre-computed.
            yn_ul_pr = MarkerGridMapping.upr_left_index_pressure(ynum, yn, y_marker, gridy_vx)
            xn_ul_pr = MarkerGridMapping.upr_left_index_pressure(xnum, xn, x_marker, gridx_vy)
            dy_ul_pr = MarkerGridMapping.upr_left_dist_pressure(yn_ul_pr, y_marker, gridy_vx, ystp_vx)
            dx_ul_pr = MarkerGridMapping.upr_left_dist_pressure(xn_ul_pr, x_marker, gridx_vy, xstp_vy)
            @inbounds begin
                # Interpolate grid normal strain rate to marker
                marker_exx[imarker] = GridToMarker.get_marker_value(yn_ul_pr, xn_ul_pr, dy_ul_pr, dx_ul_pr, exx)
                # Interpolate grid shear strain rate to marker
                marker_exy[imarker] = GridToMarker.get_marker_value(yn, xn, dy, dx, exy)
            end
        end
    end
    return nothing
end

end # module 