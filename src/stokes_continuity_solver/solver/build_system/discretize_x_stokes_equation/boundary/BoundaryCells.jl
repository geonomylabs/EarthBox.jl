module BoundaryCells

import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices
import EarthBox.Domain: outside_prescribed_velocity_zone_x_stokes

""" Calculate coefficients and rhs term for x-Stokes at boundary cells.

Boundary cells are defined as cells located along the right boundary or
cells that are inside of prescribed velocity zones.

Coefficient and right-hand side term are calculated for the large matrix
row associated with the vx unknown for the current cell with upper-left
node (i,j).

The approach described in Figure 2 of the main module docstring is used to
discretize the x-Stokes equation.

This function applies only to cells with upper left nodes located within
prescribed velocity zones or located at least one steps from the right
boundary of the basic grid.

Consider the following visual for the case where nodes are located at least
one steps from the right boundary of the basic grid:

  j=xnum-2       j=xnum-1      j=xnum
      x-------------o------------RB
                    ^
where o refers to a node for which this function is applicable, x refers
to is a non-applicable basic node and RB is a basic node located at right
boundary of the basic grid. For this case the stencil is not applied since
velocity is prescribed at the central vxC node and the matrix element at
location (ivx, ivx) is equal to 1 and the right-hand-side term is equal to
the prescribed velocity. For free slip and no slip boundary conditions vxC
is equal to zero.

As described by Gerya (2010) on page 102 the left and right parts of the
boundary condition equation vx = constant should be multiplied by a
constant Kbond to ensure that coefficients are relatively uniform and avoid
differences of many orders of magnitude:

      Kbond =  2*pscale/(xstpavr + ystpavr) =
                      4etan[1,1]/(xstpavr + ystpavr)^2

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
- `inz::Int`: 
    - Index of non-zero
      matrix arrays (Lii, Ljj, Li, Lj and Lv).
"""
function coefficients_and_rhs_terms(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    ivx = cell_indices.ivx
    
    bintern_zone = build_data.bc.bintern_zone
    bintern_velocity = build_data.bc.bintern_velocity
    pscale = build_data.pscale

    xstpavg = build_data.grid.xstpavg
    ystpavg = build_data.grid.ystpavg

    R = build_data.rhs.R

    Lv = build_data.system_vectors.Lv
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj

    Lv[inz] = 2.0 * pscale / (xstpavg + ystpavg)

    Li[inz] = ivx
    Lj[inz] = ivx

    Lii[inz] = i
    Ljj[inz] = j

    inz = inz + 1

    if outside_prescribed_velocity_zone_x_stokes(i, j, bintern_zone)
        R[ivx] = 0.0
    else
        # Internal prescribed horizontal velocity
        R[ivx] = 2.0 * pscale / (xstpavg + ystpavg) * bintern_velocity[4]
    end

    return inz
end

end # module 