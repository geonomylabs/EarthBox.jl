module MarkersToStokesArrays

import EarthBox.ModelDataContainer: ModelData
import ..BilinearNumerator
import ..BilinearAverage
import ..FlattenDensity
import ..BackupTransportArrays: backup_and_clear_stokes_arrays!

""" Interpolate Stokes-continuity marker data to grid arrays.

This function calls functions that perform the following tasks:

1. Backup and clear stokes transport arrays in preparation for
interpolation of marker properties to transport arrays. The transport 
arrays (suffix 1) must be cleared since they will be used to calculate 
weighted sums used in the denominator of the bilinear interpolation 
equations backup arrays (suffix 0) are used if no markers are found in the
surrounding grid cells (i.e. marker weights are equal to zero) during
interpolation of marker properties to transport arrays.

2. Interpolate marker information to the grid arrays for density (rho1),
viscoplastic viscosity (etan1, etas1), flow viscosity (eta_flow),
shear modulus (mus1, mun1), viscoelastic stress including rotation
(sxx1, sxy1), plastic deformation indicators (plastics, plasticn),
friction angle (fric_degrees_grid), cohesion (cohesion_grid) and
dilatation (dilatation_grid). If issues are encountered (no markers
found) then values in grids nodes are copied from backup arrays with
suffix 0 (rho0, etan0, etas0, mus0, mun0, sxx0, sxy0,
fric_degrees_grid0, cohesion_grid0, dilatation_grid0)

# Arguments
- `pymodel::PyModelData`: Main model container object.

# Returns
- `Nothing`
"""
function interpolate_markers_to_transport_arrays!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    backup_and_clear_stokes_arrays!(model)
    # Basic and pressure grids
    BilinearNumerator.calc_marker_weight_sums_for_basic_and_pressure_grids_stokes!(model, inside_flags)
    # Vy grid
    BilinearNumerator.calc_marker_weight_sums_for_vy_grid_stokes!(model, inside_flags)
    # Calculate averages on all staggered grids
    BilinearAverage.calc_marker_average_at_staggered_nodes_stokes!(model)
    # Flatten density in sticky air layer
    FlattenDensity.flatten_staggered_grid_density_in_sticky_air_layer!(model)    
    return nothing
end

""" Call functions to interpolate viscosity marker data to arrays.

# Arguments
- `model::ModelData`: Main model container object.

# Returns
- `Nothing`
"""
function interpolate_viscosity_to_transport_arrays!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    BilinearNumerator.calc_marker_weight_sums_viscosity!(model, inside_flags)
    BilinearAverage.calc_marker_average_at_nodes_viscosity!(model)
    return nothing
end

end # module 