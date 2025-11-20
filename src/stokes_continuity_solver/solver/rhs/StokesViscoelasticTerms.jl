module StokesViscoelasticTerms

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.Arrays: ArrayUtils
import ...Viscoelastic

""" Calculate viscoelastic stress and viscosity terms for Stokes equations.

Note that viscoelastic stress terms for the Stokes equations are not
viscoelastic stress. Viscoelastic stress is calculated after solving the Stokes
equation and is denoted with arrays sxx2 and sxy2.
"""

"""
    calc_viscoelastic_terms(model)

Call functions to calculate viscoelastic viscosity and stress terms.

Viscoelastic viscosity and stress term arrays (i.e. etas0, etan0, sxy0 and sxx0) 
for Stokes equations are calculated using viscosity and stress interpolated from 
markers (i.e. sxy1, sxx1, etas1, etan1) and a viscoelastic factor defined below. 
See equation 13.16 on page 185 of Gerya (2010) for more details on viscoelastic 
terms for Stokes equations.
"""
function calc_viscoelastic_terms!(model::ModelData)::Nothing
    @timeit_memit "Finished calculating viscoelastic terms" begin
        viscoelastic_shear_stress_and_viscosity_terms!(model)
        viscoelastic_normal_stress_and_viscosity_terms!(model)
    end
    return nothing
end

"""
    viscoelastic_shear_stress_and_viscosity_terms(model)

Calculate shear viscoelastic terms for Stokes equations.

# Parameters
- `etas1`: Shear viscosity interpolated to basic grid from markers in Pa.s
- `mus1`: Shear elastic modulus interpolated to basic grid from markers in Pa
- `sxy1`: Shear stress interpolated to basic grid from markers in Pa
- `viscoelastic_factor`: viscosity/(viscosity + timestep*shear_modulus)

# Updated Arrays
- `etas0`: Viscoelastic shear viscosity terms for Stokes equations in Pa.s:
  etas0 = etas1*(1 - viscoelastic_factor)
- `sxy0`: Viscoelastic shear stress terms for Stokes equations in Pa:
  sxy0 = viscoelastic_factor*sxy1
"""
function viscoelastic_shear_stress_and_viscosity_terms!(model::ModelData)::Nothing
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    etas1 = model.stokes_continuity.arrays.viscosity.etas1.array
    mus1 = model.stokes_continuity.arrays.shear_modulus.mus1.array
    sxy1 = model.stokes_continuity.arrays.stress.sxy1.array
    etas0 = model.stokes_continuity.arrays.viscosity.etas0.array
    sxy0 = model.stokes_continuity.arrays.stress.sxy0.array

    for j in 1:xnum
        for i in 1:ynum
            viscoelastic_factor = Viscoelastic.calc_viscoelastic_factor(
                i, j, etas1, mus1, timestep)
            etas0[i, j] = etas1[i, j] * (1.0 - viscoelastic_factor)
            sxy0[i, j] = sxy1[i, j] * viscoelastic_factor
        end
    end
    return nothing
end

"""
    viscoelastic_normal_stress_and_viscosity_terms(model)

Calculate normal viscoelastic terms for Stokes equations.

# Parameters
- `etan1`: Normal viscosity interpolated to pressure grid from markers in Pa.s
- `mun1`: Normal elastic modulus interpolated to pressure grid from markers in Pa
- `sxx1`: Normal stress interpolated to pressure grid from markers in Pa
- `viscoelastic_factor`: viscosity/(viscosity + timestep*shear_modulus)

# Updated Arrays
- `etan0`: Viscoelastic normal viscosity terms for Stokes equations in Pa.s:
  etan0 = etan1*(1 - viscoelastic_factor)
- `sxx0`: Viscoelastic normal stress terms for Stokes equations in Pa:
  sxx0 = viscoelastic_factor*sxx1
"""
function viscoelastic_normal_stress_and_viscosity_terms!(model::ModelData)::Nothing
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    etan1 = model.stokes_continuity.arrays.viscosity.etan1.array
    mun1 = model.stokes_continuity.arrays.shear_modulus.mun1.array
    sxx1 = model.stokes_continuity.arrays.stress.sxx1.array
    etan0 = model.stokes_continuity.arrays.viscosity.etan0.array
    sxx0 = model.stokes_continuity.arrays.stress.sxx0.array

    for j in 1:(xnum-1)
        for i in 1:(ynum-1)
            viscoelastic_factor = Viscoelastic.calc_viscoelastic_factor(
                i, j, etan1, mun1, timestep)
            etan0[i, j] = etan1[i, j] * (1 - viscoelastic_factor)
            sxx0[i, j] = sxx1[i, j] * viscoelastic_factor
        end
    end
    return nothing
end

end # module 