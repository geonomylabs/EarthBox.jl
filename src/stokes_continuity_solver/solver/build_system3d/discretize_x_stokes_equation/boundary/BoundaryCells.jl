module BoundaryCells

import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices
import ...Predicates3D: outside_prescribed_velocity_zone_x_stokes_3d

""" Calculate coefficients and rhs term for x-Stokes at boundary cells.

Boundary cells are defined as cells located along the right x-boundary
(j = xnum - 1) or cells inside a prescribed-vx zone. For these cells the
stencil is not applied since vxC is prescribed: the matrix element at
(ivx, ivx) is set to a Kbond-scaled 1 and the right-hand-side is set to
the prescribed velocity (zero for free/no slip).

The Kbond scaling described by Gerya (2010) on page 102,

    Kbond = 2*pscale / (xstpavg + ystpavg + zstpavg)

uses the average grid spacings along all three axes in 3D.

# Arguments
- `inz::Int`:
    - Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`:
    - Cell index information for the current cell.

- `build_data::StokesBuildData`:
    - Build data.

# Returns
- `inz::Int`:
    - Updated index of non-zero matrix arrays.
"""
function coefficients_and_rhs_terms(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivx = cell_indices.ivx

    bintern_zone = build_data.bc.bintern_zone
    bintern_velocity = build_data.bc.bintern_velocity
    pscale = build_data.pscale

    xstpavg = build_data.grid.xstpavg
    ystpavg = build_data.grid.ystpavg
    zstpavg = build_data.grid.zstpavg

    R = build_data.rhs.R

    Lv = build_data.system_vectors.Lv
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj

    kbond = 2.0 * pscale / (xstpavg + ystpavg + zstpavg)

    Lv[inz] = kbond
    Li[inz] = ivx
    Lj[inz] = ivx
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1

    if outside_prescribed_velocity_zone_x_stokes_3d(i, j, k, bintern_zone)
        R[ivx] = 0.0
    else
        R[ivx] = kbond * bintern_velocity[1]
    end

    return inz
end

end # module
