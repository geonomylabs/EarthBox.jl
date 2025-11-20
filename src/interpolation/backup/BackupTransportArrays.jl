module BackupTransportArrays

using EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit

"""
    backup_and_clear_stokes_arrays!(model::ModelData)

Copy the current transport arrays to backup arrays and then clear them.

This function copies Stokes-continuity transport grid arrays to backup arrays prior to clearing 
transport arrays. Clearing transport arrays is required since the arrays are used to calculate sums.

Copied grid arrays include:
- density `rho1`
- viscoplastic viscosity on basic and pressure grids `etas1` and `etan1`
- flow viscosity `eta_flow`
- shear modulus on basic and pressure grids `mus1` and `mun1`
- visco-elasto-plastic stress on basic and pressure grids `sxy1` and `sxx1`
- plastic deformation indicators on basic and pressure grids `plastics` and `plasticn`
- friction angle `fric_degrees_grid`
- cohesion `cohesion_grid`
- dilatation `dilatation_grid`

The associated backup arrays have suffix 0: `rho0`, `etan0`, `etas0`, `eta_flow0`, `mus0`, `mun0`, 
`sxx0`, `sxy0`, `plastics0`, `plasticn0`, `fric_degres_grid0`, `cohesion_grid0`, `dilatation_grid0`.

The backup transport arrays (suffix 0) are used during marker interpolation if marker weights are 
equal to zero at a grid node (i.e. no markers are present in surrounding grid cells).
"""
function backup_and_clear_stokes_arrays!(model::ModelData)
    @timeit_memit "Finished backing up and clearing Stokes transport arrays" begin
        backup_stokes_transport_arrays!(model)
        clear_stokes_transport_arrays!(model)
    end
end

"""
    backup_stokes_transport_arrays!(model::ModelData)

Copy the current transport arrays to backup arrays.

The current transport arrays are denoted with suffix "1" whereas the backup arrays are denoted with suffix "0".
"""
function backup_stokes_transport_arrays!(model::ModelData)
    model.stokes_continuity.arrays.viscosity.etas0.array = copy(model.stokes_continuity.arrays.viscosity.etas1.array)
    model.stokes_continuity.arrays.viscosity.etan0.array = copy(model.stokes_continuity.arrays.viscosity.etan1.array)
    model.stokes_continuity.arrays.viscosity.eta_flow0.array = copy(model.stokes_continuity.arrays.viscosity.eta_flow.array)
    model.stokes_continuity.arrays.shear_modulus.mus0.array = copy(model.stokes_continuity.arrays.shear_modulus.mus1.array)
    model.stokes_continuity.arrays.shear_modulus.mun0.array = copy(model.stokes_continuity.arrays.shear_modulus.mun1.array)
    model.stokes_continuity.arrays.stress.sxy0.array = copy(model.stokes_continuity.arrays.stress.sxy1.array)
    model.stokes_continuity.arrays.stress.sxx0.array = copy(model.stokes_continuity.arrays.stress.sxx1.array)
    model.stokes_continuity.arrays.density.rho0.array = copy(model.stokes_continuity.arrays.density.rho1.array)
    model.stokes_continuity.arrays.density.rho0_vy.array = copy(model.stokes_continuity.arrays.density.rho1_vy.array)
    model.stokes_continuity.arrays.plastic_def.plastics0.array = copy(model.stokes_continuity.arrays.plastic_def.plastics.array)
    model.stokes_continuity.arrays.plastic_def.plasticn0.array = copy(model.stokes_continuity.arrays.plastic_def.plasticn.array)
    model.stokes_continuity.arrays.plastic_def.cohesion_grid0.array = copy(model.stokes_continuity.arrays.plastic_def.cohesion_grid.array)
    model.stokes_continuity.arrays.plastic_def.fric_degrees_grid0.array = copy(model.stokes_continuity.arrays.plastic_def.fric_degrees_grid.array)
    model.stokes_continuity.arrays.plastic_def.dilatation_grid0.array = copy(model.stokes_continuity.arrays.plastic_def.dilatation_grid.array)
end

"""
    clear_stokes_transport_arrays!(model::ModelData)

Clear current transport arrays.

The current transport arrays are denoted with suffix "1".
"""
function clear_stokes_transport_arrays!(model::ModelData)
    fill!(model.stokes_continuity.arrays.viscosity.etas1.array, 0.0)
    fill!(model.stokes_continuity.arrays.viscosity.etan1.array, 0.0)
    fill!(model.stokes_continuity.arrays.viscosity.eta_flow.array, 0.0)
    fill!(model.stokes_continuity.arrays.shear_modulus.mus1.array, 0.0)
    fill!(model.stokes_continuity.arrays.shear_modulus.mun1.array, 0.0)
    fill!(model.stokes_continuity.arrays.stress.sxy1.array, 0.0)
    fill!(model.stokes_continuity.arrays.stress.sxx1.array, 0.0)
    fill!(model.stokes_continuity.arrays.density.rho1.array, 0.0)
    fill!(model.stokes_continuity.arrays.density.rho1_vy.array, 0.0)
    fill!(model.stokes_continuity.arrays.plastic_def.cohesion_grid.array, 0.0)
    fill!(model.stokes_continuity.arrays.plastic_def.fric_degrees_grid.array, 0.0)
    fill!(model.stokes_continuity.arrays.plastic_def.plastics.array, 0.0)
    fill!(model.stokes_continuity.arrays.plastic_def.plasticn.array, 0.0)
    fill!(model.stokes_continuity.arrays.plastic_def.dilatation_grid.array, 0.0)
end

"""
    backup_and_clear_viscosity_transport_arrays!(model::ModelData)

Copy the current transport arrays to backup arrays and then clear.

The current transport arrays, which are denoted with suffix "1", are copied to backup arrays, 
which denoted with suffix "0", before being cleared.

The backup transport arrays (suffix 0) are used during marker interpolation if marker weights 
are equal to zero at a grid node (i.e. no markers are present in surrounding grid cells.

The main transport arrays (suffix 1) are cleared prior to applying the first order bilinear 
interpolation scheme.
"""
function backup_and_clear_viscosity_transport_arrays!(model::ModelData)
    @timeit_memit "Finished clearing and backing up viscosity transport arrays" begin
        backup_viscosity_transport_arrays!(model)
        clear_viscosity_transport_arrays!(model)
    end
end

"""
    backup_viscosity_transport_arrays!(model::ModelData)

Copy the current transport arrays to backup arrays.

The current transport arrays are denoted with suffix "1" whereas the backup arrays are denoted with suffix "0".
"""
function backup_viscosity_transport_arrays!(model::ModelData)
    model.stokes_continuity.arrays.viscosity.etas0.array = copy(model.stokes_continuity.arrays.viscosity.etas1.array)
    model.stokes_continuity.arrays.viscosity.etan0.array = copy(model.stokes_continuity.arrays.viscosity.etan1.array)
end

"""
    clear_viscosity_transport_arrays!(model::ModelData)

Clear current transport arrays.

The current transport arrays are denoted with suffix "1".
"""
function clear_viscosity_transport_arrays!(model::ModelData)
    fill!(model.stokes_continuity.arrays.viscosity.etas1.array, 0.0)
    fill!(model.stokes_continuity.arrays.viscosity.etan1.array, 0.0)
end

"""
    backup_and_clear_heat_transport_arrays(model::ModelData)

Copy the current transport arrays to backup arrays and then clear.

The current transport arrays, which are denoted with suffix "1", are copied to backup arrays, 
which denoted with suffix "0", before being cleared.
"""
function backup_and_clear_heat_transport_arrays!(model::ModelData)
    backup_heat_transport_arrays!(model)
    clear_heat_transport_arrays!(model)
end

"""
    backup_heat_transport_arrays(model::ModelData)

Copy the current transport arrays to backup arrays.
"""
function backup_heat_transport_arrays!(model::ModelData)
    # Note that tk0 is copied from tk2 since tk2 is the most up-to-date
    # grid thermal solution calculated using the previous transport array as
    # the initial temperature in the heat solver
    model.heat_equation.arrays.temperature.tk0.array = copy(model.heat_equation.arrays.temperature.tk2.array)
    model.heat_equation.arrays.rhocp.rhocp0.array = copy(model.heat_equation.arrays.rhocp.rhocp1.array)
    model.heat_equation.arrays.thermal_conductivity.kt0.array = copy(model.heat_equation.arrays.thermal_conductivity.kt1.array)
    model.heat_equation.arrays.radiogenic_production.hr0.array = copy(model.heat_equation.arrays.radiogenic_production.hr1.array)
    model.heat_equation.arrays.adiabatic_production.ha0.array = copy(model.heat_equation.arrays.adiabatic_production.ha1.array)
end

"""
    clear_heat_transport_arrays!(model::ModelData)

Clear current transport arrays.

The current transport arrays are denoted with suffix "1".
"""
function clear_heat_transport_arrays!(model::ModelData)
    fill!(model.heat_equation.arrays.temperature.tk1.array, 0.0)
    fill!(model.heat_equation.arrays.rhocp.rhocp1.array, 0.0)
    fill!(model.heat_equation.arrays.thermal_conductivity.kt1.array, 0.0)
    fill!(model.heat_equation.arrays.radiogenic_production.hr1.array, 0.0)
    fill!(model.heat_equation.arrays.adiabatic_production.ha1.array, 0.0)
end

end # module 