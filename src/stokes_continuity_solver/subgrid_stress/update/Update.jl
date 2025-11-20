module Update

import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.Interpolation: BilinearAverage
import EarthBox.ModelDataContainer: ModelData
import ..SubgridShearStress
import ..SubgridNormalStress
import ..RemainingStress

""" Calculate subgrid stress change and apply correction to marker stress.

# Algorithm
1. Interpolate pre-Stokes-solver stress (sxy1, sxx1) back to markers. These
   interpolated nodal stresses are sxy_marker_from_grid and sxx_marker_from_grid.

2. Calculate the relaxed subgrid stress (sxy_marker_diff_relax, 
   sxx_marker_diff_relax) by taking the difference between the interpolated
   nodal grid stress (sxy_marker_from_grid, sxx_marker_from_grid) and the
   current marker stress (marker_sxy, marker_sxx):

   sxy_marker_diff_relax =
       (sxy_marker_from_grid - marker_sxy)*stress_diffusion_term,
   sxx_marker_diff_relax =
       (sxx_marker_from_grid - marker_sxx)*stress_diffusion_term,

   The stress_diffusion_term is defined as follows:

   stress diffusion_term =
       (1 - exp(-subgrid_diff_coef_stress*timestep/maxwell_time))

   where subgrid_diff_coef_stress is the diffusion coefficient, the
   maxwell_time = marker_eta/shear_modulus,
   where marker_eta is the viscoplastic viscosity of the marker and
   shear_modulus is the shear modulus of the marker.

3. Add the relaxed subgrid marker stress to the current marker
   stress arrays (marker_sxy, marker_sxx):
       marker_sxy = marker_sxy + sxy_marker_diff_relax,
       marker_sxx = marker_sxx + sxx_marker_diff_relax,

4. Calculate the relaxed subgrid stress on grid (dsxyn, dsxxn)
   by interpolating relaxed subgrid marker stress
   (sxy_marker_diff_relax, sxx_marker_diff_relax) from markers to grid
   arrays.

5. Calculate the remaining stress change (dsxy, dsxx) by removing the
   the relaxed subgrid stress (dsxyn, dsxxn) from the viscoelastic stress 
   change associated with the updated Stokes solution (dsxy, dsxx):
       dsxy = dsxy - dsxyn
       dsxx = dsxx - dsxxn

# Background
Stress changes are decomposed into a subgrid part and a remaining part:

dstress_grid = dstress_subgrid_grid + dstress_remaining_grid

The subgrid stress is calculated by applying a time-dependent relaxation
term by the difference between stress interpolated from the grid to
markers and current marker stress:

dstress_subgrid_marker =
        (stress_interp_from_grid - marker_stress)*relaxation_term

The remaining part is calculated subtracting the subgrid stress from the
viscoelastic stress computing using the updated Stokes solution:

dstress_remaining_grid = dstress_viscoelastic_grid - dstress_subgrid_grid
"""
function update_subgrid_stress!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "Finished updating subgrid stress" begin
        # Initialize arrays to zero
        model.stokes_continuity.arrays.stress_change.dsxyn.array .= 0.0
        model.stokes_continuity.arrays.stress_change.dsxxn.array .= 0.0
        model.interpolation.arrays.grid_weights.wtetas.array .= 0.0
        model.interpolation.arrays.grid_weights.wtetan.array .= 0.0
        # Calculate subgrid stress changes
        SubgridShearStress.subgrid_shear_stress_changes!(model, inside_flags)
        SubgridNormalStress.subgrid_normal_stress_changes!(model, inside_flags)
        BilinearAverage.calc_marker_average_subgrid_stress_change!(model)
        # Calculate remaining stress change
        RemainingStress.calculate_remaining_stress_change!(model)
    end
    return nothing
end

end # module 