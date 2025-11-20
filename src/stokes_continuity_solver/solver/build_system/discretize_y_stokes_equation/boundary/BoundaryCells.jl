module BoundaryCells

import EarthBox.Domain: outside_prescribed_velocity_zone_y_stokes
import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices

""" Calculate coefficients and rhs term for y-Stokes at boundary cells.

Boundary cells are defined as cells located along the lower boundary or
cells that are inside of prescribed velocity zones.

Coefficient and right-hand side term are calculated for the large matrix
row associated with the vy unknown for the current cell with upper-left
node (i,j).

The approach described in Figure 3 of the main module docstring is used to
discretize the y-Stokes equation.

The function applies only to basic cells with upper left node located at
least one step from the bottom boundary of the basic grid. Consider the
following visual:

      x i=ynum-2
      |
      |
      |
      o i=ynum-1
      |
      |
      |
     LB i=ynum

where o refers to a node for which this function is applicable, x refers to
a non-applicable basic node and LB is a basic node located at lower
boundary of the basic grid.

For this case the stencil is not applied since velocity is prescribed at
the central vyC node and the matrix element at location (ivy, ivy) is equal
to 1 and the right-hand-side term is equal to the prescribed velocity. For
free slip and no slip boundary conditions vyC is equal to zero.

As described by Gerya (2010) on page 102 the left and right parts of the
boundary condition equation vy = constant should be multiplied by a
constant Kbond to ensure that coefficients are relatively uniform and avoid
differences of many orders of magnitude:

      Kbond =  2*pscale/(xstpavr+ystpavr) =
                      4etan[1,1]/(xstpavr+ystpavr)^2

This coefficient scaling facilitates a more accurate solution when applying
direct solvers.

# Arguments
- `inz::Int`: 
    - Index of non-zero
      matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`: 
    - Cell index information for the current cell.

- `build_data::StokesBuildData`: 
    - Build data.

# Returns
- `inz::Int`: Index of non-zero
  matrix arrays (Lii, Ljj, Li, Lj and Lv).
"""
function coefficients_and_rhs_terms(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    ivy = cell_indices.ivy
    xstpavg = build_data.grid.xstpavg
    ystpavg = build_data.grid.ystpavg
    pscale = build_data.pscale
    bintern_zone = build_data.bc.bintern_zone
    bintern_velocity = build_data.bc.bintern_velocity
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    Lv[inz] = 2.0 * pscale / (xstpavg + ystpavg)
    #Lv[inz]=1.0 # in Stokes_Continuity_solver_channel.m this is set to 1
    Li[inz] = ivy
    Lj[inz] = ivy
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1
    if outside_prescribed_velocity_zone_y_stokes(i, j, bintern_zone)
        R[ivy] = 0.0
    else
        # Internal prescribed horizontal velocity
        R[ivy] = 2.0 * pscale / (xstpavg + ystpavg) * bintern_velocity[8]
    end
    return inz
end

end # module 