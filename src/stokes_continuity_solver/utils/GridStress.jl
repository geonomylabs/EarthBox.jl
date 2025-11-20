module GridStress

import EarthBox.ModelDataContainer: ModelData
import ..Viscoelastic

"""
    forecast_viscoelastic_grid_stress!(model::ModelData)

Forecast viscoelastic shear (sxy2) and normal stress (sxx2) on grid.

Viscoelastic shear stress on the basic grid is calculated using the following 
equation:

    sxy2[i,j] = (1.0 - viscoelastic_factor_xy)*2.0*etas1[i,j]*exy[i,j] + 
                 viscoelastic_factor_xy*sxy1[i,j]

where etas1 and sxy1 are the current shear viscoplastic viscosity and shear stress,
respectively, exy is the shear strain rate calculated using the updated Stokes 
solution, and viscoelastic_factor_xy is given by:

    viscoelastic_factor_xy = etas1[i,j]/(etas1[i,j] + timestep*mus1[i,j])

where mus1 is the elastic shear modulus and timestep is the main model timestep.

Viscoelastic normal stress sxx1 on the pressure grid is calculated using a similar 
approach but with normal viscosity and normal stress interpolated from markers 
(etan1, sxx1), the normal strain rate from the updated Stokes solution (exx), and 
the normal elastic modulus interpolated from markers (mun1).
"""
function forecast_viscoelastic_grid_stress!(model::ModelData)
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    forecast_viscoelastic_shear_stress!(model, timestep, ynum, xnum)
    forecast_viscoelastic_normal_stress!(model, timestep, ynum, xnum)
end

"""
    forecast_viscoelastic_shear_stress!(model::ModelData, timestep::Float64, 
                                     ynum::Int, xnum::Int)

Calculate viscoelastic shear stress on grid nodes.

# Arguments
- `timestep`: Main model time step in seconds
- `ynum`: Number of grid basic nodes in y-direction
- `xnum`: Number of basic grid nodes in x-direction

# Updated Arrays
- `model.stokes_continuity.arrays.stress.sxy2`: Viscoelastic shear stress (Pa) on 
  basic grid
"""
function forecast_viscoelastic_shear_stress!(
    model::ModelData, timestep::Float64, ynum::Int, xnum::Int
)::Nothing
    etas1 = model.stokes_continuity.arrays.viscosity.etas1.array
    mus1 = model.stokes_continuity.arrays.shear_modulus.mus1.array
    exy = model.stokes_continuity.arrays.strain_rate_and_spin.exy.array
    sxy1 = model.stokes_continuity.arrays.stress.sxy1.array
    sxy2 = model.stokes_continuity.arrays.stress.sxy2.array

    fill!(sxy2, 0.0)
    for j in 1:xnum
        for i in 1:ynum
            viscoelastic_factor = Viscoelastic.calc_viscoelastic_factor(i, j, etas1, mus1, timestep)
            sxy2[i,j] =(
                (1.0 - viscoelastic_factor) * 2.0 * etas1[i,j] * exy[i,j] 
                + viscoelastic_factor*sxy1[i,j]
            )
        end
    end
    return nothing
end

"""
    forecast_viscoelastic_normal_stress!(model::ModelData, timestep::Float64, 
                                      ynum::Int, xnum::Int)

Calculate viscoelastic normal stress on grid nodes.

# Arguments
- `timestep`: Main model time step in seconds
- `ynum`: Number of grid basic nodes in y-direction
- `xnum`: Number of basic grid nodes in x-direction

# Updated Arrays
- `model.stokes_continuity.arrays.stress.sxx2`: Viscoelastic normal stress (Pa) on 
  pressure grid
"""
function forecast_viscoelastic_normal_stress!(
    model::ModelData, timestep::Float64, ynum::Int, xnum::Int
)::Nothing
    etan1 = model.stokes_continuity.arrays.viscosity.etan1.array
    mun1 = model.stokes_continuity.arrays.shear_modulus.mun1.array
    exx = model.stokes_continuity.arrays.strain_rate_and_spin.exx.array
    sxx1 = model.stokes_continuity.arrays.stress.sxx1.array
    sxx2 = model.stokes_continuity.arrays.stress.sxx2.array

    fill!(sxx2, 0.0)
    for j in 1:(xnum-1)
        for i in 1:(ynum-1)
            viscoelastic_factor = Viscoelastic.calc_viscoelastic_factor(i, j, etan1, mun1, timestep)
            sxx2[i,j] = (
                (1.0 - viscoelastic_factor) * 2.0 * etan1[i,j] * exx[i,j] 
                + viscoelastic_factor*sxx1[i,j]
                )
        end
    end
    return nothing
end

"""
    calculate_deviatoric_grid_stress_change!(model::ModelData)

Calculate grid normal and shear stress changes.

Normal (dsxx) and shear (dsxy) stress changes are calculated using old stress 
interpolated from markers (sxx1, sxy1) and new viscoelastic stress (sxx2, sxy2) 
based on updated strain rates from new Stokes velocity solution.

# Updated Arrays
- `model.stokes_continuity.arrays.stress_change.dsxy`: Shear stress change on 
  basic grid
- `model.stokes_continuity.arrays.stress_change.dsxx`: Normal stress change on 
  pressure grid
"""
function calculate_deviatoric_grid_stress_change!(model::ModelData)
    stress = model.stokes_continuity.arrays.stress
    stress_change = model.stokes_continuity.arrays.stress_change
    dsxy = stress.sxy2.array - stress.sxy1.array
    dsxx = stress.sxx2.array - stress.sxx1.array
    stress_change.dsxy.array = copy(dsxy)
    stress_change.dsxx.array = copy(dsxx)
end

end # module 