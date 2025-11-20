module SubgridShearStress

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Interpolation: GridToMarker
import EarthBox.Interpolation: WeightFuncs
import EarthBox: GridFuncs
import ..DiffusionTerm

""" Calculate subgrid shear stress changes and correct marker stress.

# Algorithm
1. Interpolate to markers the stress from basic grid (i.e. sxy1). Note that sxy1 was 
   interpolated from markers prior to Stokes solution.
2. Calculate the nodal-marker stress difference between the interpolated grid stress and 
   the current marker stress and apply subgrid relaxation to nodal-marker stress difference.
3. Add the relaxed stress difference to the current marker stress.
4. Update grid numerator and denominator sums for interpolation of relaxed subgrid marker 
   stress difference to grid.

# Updated Arrays
- model.markers.arrays.stress.marker_sxy: Marker shear stress
- model.stokes_continuity.arrays.stress_change.dsxyn: Subgrid nodal-marker shear stress 
  change (Pa) on basic grid
- model.interpolation.arrays.grid_weights.wtetas: Marker interpolation weight sum on 
  basic grid
"""
function subgrid_shear_stress_changes!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    subgrid_diff_coef_stress = model.markers.parameters.subgrid_diffusion.subgrid_diff_coef_stress.value
    
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_xn = model.markers.arrays.grid_marker_relationship.marker_xn.array
    marker_yn = model.markers.arrays.grid_marker_relationship.marker_yn.array
    marker_dx = model.markers.arrays.grid_marker_relationship.marker_dx.array
    marker_dy = model.markers.arrays.grid_marker_relationship.marker_dy.array
    marker_wtforULnode = model.interpolation.arrays.marker_weights.marker_wtforULnode.array
    marker_wtforLLnode = model.interpolation.arrays.marker_weights.marker_wtforLLnode.array
    marker_wtforURnode = model.interpolation.arrays.marker_weights.marker_wtforURnode.array
    marker_wtforLRnode = model.interpolation.arrays.marker_weights.marker_wtforLRnode.array
    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_sxy = model.markers.arrays.stress.marker_sxy.array
    mat_mu = model.materials.arrays.mat_mu.array
    sxy1 = model.stokes_continuity.arrays.stress.sxy1.array
    dsxyn = model.stokes_continuity.arrays.stress_change.dsxyn.array
    wtetas = model.interpolation.arrays.grid_weights.wtetas.array

    for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                material_id = marker_matid[imarker]
                viscosity = marker_eta[imarker]
                shear_modulus = mat_mu[material_id]
                # Define interpolation parameters
                dx_upr_left = marker_dx[imarker]
                dy_upr_left = marker_dy[imarker]
                ix_upr_left = marker_xn[imarker]
                iy_upr_left = marker_yn[imarker]
                marker_weight_upr_left = marker_wtforULnode[imarker]
                marker_weight_lower_left = marker_wtforLLnode[imarker]
                marker_weight_upr_right = marker_wtforURnode[imarker]
                marker_weight_lower_right = marker_wtforLRnode[imarker]
            end
            stress_diffusion_term = DiffusionTerm.calculate_subgrid_stress_diffusion_term(
                viscosity, shear_modulus, timestep, subgrid_diff_coef_stress)
            sxy_marker_from_grid = GridToMarker.get_marker_value_from_basic_grid_array(
                iy_upr_left, ix_upr_left, dy_upr_left, dx_upr_left, sxy1
            )
            # Calculate nodal-marker relaxed subgrid stress difference
            sxy_marker_diff_relax = (
                (sxy_marker_from_grid - marker_sxy[imarker]) * stress_diffusion_term
            )
            # Correcting stress for the marker
            @inbounds marker_sxy[imarker] += sxy_marker_diff_relax
            # Update numerator sum for interpolation of relaxed marker stress difference 
            # back to grid (dsxyn)
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                dsxyn, sxy_marker_diff_relax
            )
            WeightFuncs.update_weight_exclusive!(
                iy_upr_left, ix_upr_left,
                marker_weight_upr_left, marker_weight_lower_left,
                marker_weight_upr_right, marker_weight_lower_right, 1.0,
                dx_upr_left, dy_upr_left,
                wtetas, 1.0
            )
        end
    end
end

end # module 