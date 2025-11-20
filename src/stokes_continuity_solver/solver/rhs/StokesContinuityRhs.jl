module StokesContinuityRhs

using Printf
import EarthBox: Domain
import EarthBox.Arrays: ArrayUtils
import EarthBox.PrintFuncs: @timeit_memit

"""
    viscoelastic_rhs!(model)

Call functions to calculate rhs components of Stokes-continuity equations.

The right-hand side terms for the Stokes equations are a function of gravitational 
acceleration (gravity_y and gravity_x), density interpolated from markers (rho1) and 
stress multiplied by a visco-elasticity factor (sxx0 and sxy0).

We refer to the product of stress and the visco-elasticity factor as a viscoelastic 
stress term. Viscoelastic stress terms sxx0 and sxy0 are calculated in a previous step. 
See equation 13.16 on page 185 of Gerya (2010) for details on right-hand side terms in 
the Stokes equations including viscoelastic stress terms.

# Updated Array Objects
- `RX1`: Array((ynum+1, xnum), Float64) - Right-hand part for x-Stokes equation on 
  staggered vx grid.
- `RY1`: Array((ynum, xnum+1), Float64) - Right-hand part for y-Stokes equation on 
  staggered vy grid.
- `RC1`: Array((ynum-1, xnum-1), Float64) - Right-hand part for continuity equation on 
  pressure grid (currently set to zero).
"""
function viscoelastic_rhs!(model)
    @timeit_memit "Finished calculating rhs components for Stokes-continuity" begin
      ArrayUtils.setzeros!(model.stokes_continuity.arrays.rhs.RX1)
      ArrayUtils.setzeros!(model.stokes_continuity.arrays.rhs.RY1)
      ArrayUtils.setzeros!(model.stokes_continuity.arrays.rhs.RC1)

      viscoelastic_rhs_y_stokes!(model)
      viscoelastic_rhs_x_stokes!(model)
      viscoelastic_rhs_continuity!(model)
    end
end

"""
    viscoelastic_rhs_y_stokes!(model)

Compute right-hand-side part of y-Stokes-continuity equations.

# Updated Array Objects
- `RY1`: Array((ynum, xnum+1), Float64) - Right-hand part for y-Stokes equation on 
  staggered vy grid.

# Parameters
- `xnum`, `ynum`: Int, Int - Number of grid nodes in x and y-directions respectively.
- `gravity_y`: Float64 - Y-component of gravitational acceleration in m/s/s.
- `xstp_b`: Array((xnum-1), Float64) - Width of cells in x-direction in meters for basic grid.
- `ystp_vx`: Array((ynum), Float64) - Width of staggered Vx cells in y-direction in meters.
- `rho1`: Array((ynum, xnum), Float64) - Density interpolated from markers in kg/m^3 
  defined on basic grid.
- `sxx0`: Array((ynum-1, xnum-1), Float64) - Viscoelastic normal stress terms for Stokes 
  equations in Pa defined on pressure grid.
- `sxy0`: Array((ynum, xnum), Float64) - Viscoelastic shear stress terms for Stokes 
  equations in Pa defined on basic grid.
"""
function viscoelastic_rhs_y_stokes!(model)
    gravity_y = model.gravity.parameters.gravity_y.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    xstp_b = model.grids.arrays.basic.xstp_b.array
    ystp_vx = model.grids.arrays.staggered_vx.ystp_vx.array
    rho1 = model.stokes_continuity.arrays.density.rho1.array
    sxx0 = model.stokes_continuity.arrays.stress.sxx0.array
    sxy0 = model.stokes_continuity.arrays.stress.sxy0.array
    RY1 = model.stokes_continuity.arrays.rhs.RY1.array

    # Density array defined on vy grid used for drunken sailor stabilization
    rho1_vy = model.stokes_continuity.arrays.density.rho1_vy.array
    _iuse_interface_stabilization = model.stokes_continuity.parameters.build.
        iuse_interface_stabilization.value

    # Looping over staggered Vy grid nodes
    for j in 2:xnum
        for i in 2:ynum
            if !Domain.basic_node_along_lower_boundary(i, ynum)
                # Adding gravity term
                if _iuse_interface_stabilization == 1
                    # Use density defined on vy grid for drunken sailor interface 
                    # stabilization
                    RY1[i,j] = -gravity_y * rho1_vy[i,j]
                else
                    # Use density defined on basic grid (no interface stabilization)
                    RY1[i,j] = -gravity_y * (rho1[i,j] + rho1[i,j-1]) / 2
                end
                # Adding dsyy0/dy where syy0 = -sxx0 (deviatoric stress)
                RY1[i,j] += (sxx0[i,j-1] - sxx0[i-1,j-1]) / ystp_vx[i]
                # Adding dsyx0/dx where syx0 = sxy0
                RY1[i,j] -= (sxy0[i,j] - sxy0[i,j-1]) / xstp_b[j-1]
            end
        end
    end
end

"""
    viscoelastic_rhs_x_stokes(model)

Compute right-hand-side part of x-Stokes-continuity equations.

# Updated Array Objects
- `RX1`: Array((ynum+1, xnum), Float64) - Right-hand part for x-Stokes equation on 
  staggered vx grid.

# Parameters
- `xnum`, `ynum`: Int, Int - Number of grid nodes in x and y-directions respectively.
- `gravity_x`: Float64 - X-component of gravitational acceleration in m/s/s.
- `ystp_b`: Array((ynum-1), Float64) - Width of cells in y-direction in meters for basic grid.
- `xstp_vy`: Array((xnum), Float64) - Width of staggered Vy cells in x-direction in meters.
- `rho1`: Array((ynum, xnum), Float64) - Density interpolated from markers in kg/m^3 
  defined on basic grid.
- `sxx0`: Array((ynum-1, xnum-1), Float64) - Viscoelastic normal stress terms for Stokes 
  equations in Pa defined on pressure grid.
- `sxy0`: Array((ynum, xnum), Float64) - Viscoelastic shear stress terms for Stokes 
  equations in Pa defined on basic grid.
"""
function viscoelastic_rhs_x_stokes!(model)
    gravity_x = model.gravity.parameters.gravity_x.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    ystp_b = model.grids.arrays.basic.ystp_b.array
    xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array
    rho1 = model.stokes_continuity.arrays.density.rho1.array
    sxx0 = model.stokes_continuity.arrays.stress.sxx0.array
    sxy0 = model.stokes_continuity.arrays.stress.sxy0.array
    RX1 = model.stokes_continuity.arrays.rhs.RX1.array

    for j in 2:xnum
    	for i in 2:ynum
            if !Domain.basic_node_on_right_boundary(j, xnum)
                # Adding gravity term
                RX1[i,j] = -gravity_x * (rho1[i,j] + rho1[i-1,j]) / 2
                # Adding dsxx0/dx
                RX1[i,j] -= (sxx0[i-1,j] - sxx0[i-1,j-1]) / xstp_vy[j]
                # Adding dsxy0/dy
                RX1[i,j] -= (sxy0[i,j] - sxy0[i-1,j]) / ystp_b[i-1]
            end
        end
    end
end

"""
    viscoelastic_rhs_continuity(model)

Compute right-hand-side part of continuity equation.

# Updated Array Objects
- `RC1`: Array((ynum-1, xnum-1), Float64) - Right-hand part for continuity equation on 
  pressure grid.
"""
function viscoelastic_rhs_continuity!(model)
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    strain_rate_and_spin = model.stokes_continuity.arrays.strain_rate_and_spin
    eii_plastic_pressure = strain_rate_and_spin.eii_plastic_pressure.array
    RC1 = model.stokes_continuity.arrays.rhs.RC1.array

    dilatation_grid = model.stokes_continuity.arrays.plastic_def.dilatation_grid.array

    for j in 1:xnum-1
        for i in 1:ynum-1
            dilatation_angle = dilatation_grid[i,j]
            plastic_strain_rate = eii_plastic_pressure[i,j]
            plastic_volumetric_effect = calc_plastic_volumetric_effect(
                dilatation_angle, plastic_strain_rate)
            RC1[i,j] = plastic_volumetric_effect
        end
    end
end

"""
    calc_plastic_volumetric_effect(dilatation_angle, plastic_strain_rate)

Calculate the volumetric effect of plastic strain.
"""
function calc_plastic_volumetric_effect(
    dilatation_angle::Float64,
    plastic_strain_rate::Float64
)::Float64
    dilatation_angle_radians = dilatation_angle * π / 180.0
    return 2.0 * sin(dilatation_angle_radians) * plastic_strain_rate
end

end # module 