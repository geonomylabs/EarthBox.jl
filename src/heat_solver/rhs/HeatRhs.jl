module HeatRhs

import EarthBox.ModelDataContainer: ModelData

function initialize_rhs_terms_array_for_heat!(model::ModelData)
    model.heat_equation.arrays.rhs.RT1.array .= 
        model.heat_equation.arrays.radiogenic_production.hr1.array
end

""" Initialize rhs with adiabatic and friction components.

Adiabatic heating is calculated using the following equation:
    alpha*T*DP/dt

where alpha is thermal expansivity (1/K), T is temperature (K) and

DP/dt ~ vx*gravity_x*rho + vy*gravity_y*rho

where vx and vy are velocity components (m/s), rho is density (kg/m3) and
gravity_x and gravity_y are gravitational acceleration components (m/s/s).

Note that the array ha1 stores alpha*T terms on grid nodes.

# Updated Arrays
- `model.heat_equation.arrays.rhs.RT1.array`: Right-hand side array terms for heat equation 
  defined on basic grid.
"""
function add_adiabatic_terms_to_rhs_grid!(model::ModelData)
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    gravity_x = model.gravity.parameters.gravity_x.value
    gravity_y = model.gravity.parameters.gravity_y.value
    vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array
    ha1 = model.heat_equation.arrays.adiabatic_production.ha1.array
    rho1 = model.stokes_continuity.arrays.density.rho1.array
    RT1 = model.heat_equation.arrays.rhs.RT1.array

    for j in 2:xnum-1
        for i in 2:ynum-1
            RT1[i, j] += ha1[i, j] * rho1[i, j] * (
                gravity_x * (vx1[i, j] + vx1[i+1, j]) +
                gravity_y * (vy1[i, j] + vy1[i, j+1])
            ) / 2.0
        end
    end
end

""" Initialize rhs grid with friction components.

Viscoelastic shear heating for Temperature nodes is calculated using the following equation:

    Hs = 2*Sxx*Sxx/2/etan + 2*Sxy*Sxy/2/etas

# Updated Arrays
- `model.heat_equation.arrays.rhs.RT1.array`: Right-hand side array terms for heat equation 
  defined on basic grid.
"""
function add_friction_terms_to_rhs_grid!(model::ModelData)
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    sxy2 = model.stokes_continuity.arrays.stress.sxy2.array
    etas1 = model.stokes_continuity.arrays.viscosity.etas1.array
    sxx2 = model.stokes_continuity.arrays.stress.sxx2.array
    etan1 = model.stokes_continuity.arrays.viscosity.etan1.array
    RT1 = model.heat_equation.arrays.rhs.RT1.array

    for j in 2:xnum-1
        for i in 2:ynum-1
            RT1[i, j] += sxy2[i, j]^2.0 / etas1[i, j]
            RT1[i, j] += (
                sxx2[i-1, j-1]^2.0 / etan1[i-1, j-1] +
                sxx2[i, j-1]^2.0 / etan1[i, j-1] +
                sxx2[i-1, j]^2.0 / etan1[i-1, j] +
                sxx2[i, j]^2.0 / etan1[i, j]
            ) / 4.0
        end
    end
end

end # module 