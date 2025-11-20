module SubgridNormalStress

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Interpolation: GridToMarker
import EarthBox.Interpolation: WeightFuncs
import EarthBox: GridFuncs
import ..DiffusionTerm

""" Calculate subgrid normal stress changes and correct marker stress.

# Algorithm
1. Interpolate to markers the stress from basic grid (i.e. sxx1). Note that sxx1 was 
   interpolated from markers prior to Stokes solution.
2. Calculate the nodal-marker stress difference between the interpolated grid stress and 
   the current marker stress and apply subgrid relaxation to nodal-marker stress difference.
3. Add the relaxed stress difference to the current marker stress.
4. Update grid numerator and denominator sums for interpolation of relaxed subgrid marker 
   stress difference to grid.

# Updated Arrays
- model.markers.arrays.stress.marker_sxx: Marker normal stress
- model.stokes_continuity.arrays.stress_change.dsxxn: Subgrid nodal-marker normal stress 
  change (Pa) on pressure grid
- model.interpolation.arrays.grid_weights.wtetan: Marker interpolation weight sum on 
  pressure grid
"""
function subgrid_normal_stress_changes!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    subgrid_diff_coef_stress = model.markers.parameters.subgrid_diffusion.subgrid_diff_coef_stress.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_wtforCnode = model.interpolation.arrays.marker_weights.marker_wtforCnode.array
    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_sxx = model.markers.arrays.stress.marker_sxx.array
    mat_mu = model.materials.arrays.mat_mu.array
    gridx_vy = model.grids.arrays.staggered_vy.gridx_vy.array
    xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array
    gridy_vx = model.grids.arrays.staggered_vx.gridy_vx.array
    ystp_vx = model.grids.arrays.staggered_vx.ystp_vx.array
    sxx1 = model.stokes_continuity.arrays.stress.sxx1.array
    dsxxn = model.stokes_continuity.arrays.stress_change.dsxxn.array
    wtetan = model.interpolation.arrays.grid_weights.wtetan.array

    for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                x_marker = marker_x[imarker]
                y_marker = marker_y[imarker]
                material_id = marker_matid[imarker]
                viscosity = marker_eta[imarker]
                sxx_marker = marker_sxx[imarker]
                shear_modulus = mat_mu[material_id]
                ix_upr_left = marker_xn[imarker]
                iy_upr_left = marker_yn[imarker]
                marker_weight_central = marker_wtforCnode[imarker]
            end
            stress_diffusion_term = DiffusionTerm.calculate_subgrid_stress_diffusion_term(
                viscosity, shear_modulus, timestep, subgrid_diff_coef_stress)
            sxx_marker_from_grid = GridToMarker.get_marker_value_from_pressure_grid_array(
                ynum, xnum, iy_upr_left, ix_upr_left, y_marker, x_marker,
                gridy_vx, ystp_vx, gridx_vy, xstp_vy, sxx1
            )
            # Calculate nodal-marker subgrid sxx stress difference
            sxx_marker_diff_relax = (sxx_marker_from_grid - sxx_marker) * stress_diffusion_term
            # Correcting old stress for the marker
            @inbounds marker_sxx[imarker] += sxx_marker_diff_relax
            # Update numerator sum for interpolation of relaxed marker stress difference 
            # back to grid (dsxxn)
            WeightFuncs.update_weight_central!(
                iy_upr_left, ix_upr_left, marker_weight_central, 1.0,
                dsxxn, sxx_marker_diff_relax
            )
            WeightFuncs.update_weight_central!(
                iy_upr_left, ix_upr_left, marker_weight_central, 1.0,
                wtetan, 1.0
            )
        end
    end
end

end # module 