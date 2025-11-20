module MarkerToHeatArrays

import EarthBox.ModelDataContainer: ModelData
import ..BilinearNumerator
import ..BilinearAverage
import ..BackupTransportArrays
import ..ApplyThermalBC

""" Interpolate marker information to heat equation transport grid arrays.

This function manages the interpolation of marker temperature, conductivity,
heat capacity, heat production, and adiabatic heat production to transport 
arrays on the basic grid. Transport grid arrays for the heat 
equation include temperature `tk1`, thermal conductivity `kt1`, the 
product of density and heat capacity `rhocp1`, heat production `hr1` and
the product thermal expansivity and temperature `ha1`.

Prior to interpolation transport arrays are backup up and cleared. Backup 
arrays have suffix 0: temperature `tk0`, thermal conductivity `kt0`, the
product of density and heat capacity `rhocp0`, heat production `hr0` and
the product thermal expansivity and temperature `ha0`. Backup arrays are
used to define values in grid nodes if there are no markers in the
corresponding grid cell during interpolation.

Boundary conditions for the heat equation are applied after
interpolation to ensure that boundary condition values are correct along
boundary nodes.
"""
function interpolate_markers_to_transport_arrays!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    BackupTransportArrays.backup_and_clear_heat_transport_arrays!(model)
    BilinearNumerator.calc_marker_weight_sums_for_basic_grid_heat!(model, inside_flags)
    BilinearAverage.calc_marker_average_at_basic_nodes_heat!(model)
    ApplyThermalBC.apply_thermal_bc_along_boundaries!(model)
    return nothing
end

end # module 