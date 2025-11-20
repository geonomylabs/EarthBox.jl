module StressCorrection

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs
import EarthBox.Interpolation: GridToMarker

"""
    apply_marker_stress_correction_for_remaining_grid_stress!(model)

Correct marker stress for remaining stress.


The total marker stress correction is defined as follows:

    marker_stress_corrected =
            marker_stress + dstress_marker_subgrid + dtress_marker_remaining

Note that the subgrid stress correction was applied in a previous step.
The current function only corrects for the remaining grid stress change:

    marker_sxy[imarker] = marker_sxy[imarker] + dsxy_remaining
    marker_sxx[imarker] = marker_sxx[imarker] + dsxx_remaining

Note that remaining grid stress change arrays dsxy and dsxx are defined as follows:
    dstress_remaining =
          dstress_viscoelastic_current - dstress_viscoelastic_old
        - dtress_subgrid

# Updated Arrays

## model.markers.arrays.stress
- marker_sxy.array: Array{Float64,1} - Marker shear stress
- marker_sxx: Array{Float64,1} - Marker normal stress
"""
function apply_marker_stress_correction_for_remaining_grid_stress!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array
    gridx_vy = model.grids.arrays.staggered_vy.gridx_vy.array
    xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array
    gridy_vx = model.grids.arrays.staggered_vx.gridy_vx.array
    ystp_vx = model.grids.arrays.staggered_vx.ystp_vx.array
    dsxy = model.stokes_continuity.arrays.stress_change.dsxy.array
    dsxx = model.stokes_continuity.arrays.stress_change.dsxx.array
    marker_sxy = model.markers.arrays.stress.marker_sxy.array
    marker_sxx = model.markers.arrays.stress.marker_sxx.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                x_marker = marker_x[imarker]
                y_marker = marker_y[imarker]
                sxy_marker = marker_sxy[imarker]
                sxx_marker = marker_sxx[imarker]
                dx_upr_left = marker_dx[imarker]
                dy_upr_left = marker_dy[imarker]
                ix_upr_left = marker_xn[imarker]
                iy_upr_left = marker_yn[imarker]
            end
            dsxy_marker = GridToMarker.get_marker_value_from_basic_grid_array(
                iy_upr_left, ix_upr_left, dy_upr_left, dx_upr_left, dsxy
            )
            @inbounds marker_sxy[imarker] = sxy_marker + dsxy_marker
            dsxx_marker = GridToMarker.get_marker_value_from_pressure_grid_array(
                ynum, xnum, iy_upr_left, ix_upr_left, y_marker, x_marker,
                gridy_vx, ystp_vx, gridx_vy, xstp_vy, dsxx
            )
            @inbounds marker_sxx[imarker] = sxx_marker + dsxx_marker
        end
    end
    return nothing
end

end # module 