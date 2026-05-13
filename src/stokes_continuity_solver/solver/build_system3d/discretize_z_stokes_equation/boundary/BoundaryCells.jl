module BoundaryCells

import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices
import ...Predicates3D: outside_prescribed_velocity_zone_z_stokes_3d

""" Calculate coefficients and rhs term for z-Stokes at boundary cells.

Boundary cells are defined as cells located along the back z-boundary
(k = znum - 1) or cells inside a prescribed-vz zone. For these cells the
stencil is not applied since vzC is prescribed: the matrix element at
(ivz, ivz) is set to a Kbond-scaled 1 and the right-hand-side is set to
the prescribed velocity (zero for free/no slip).

# Arguments
- `inz::Int`:
    - Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`:
    - Cell index information for the current cell.

- `build_data::StokesBuildData`:
    - Build data.

# Returns
- `inz::Int`: Updated index of non-zero matrix arrays.
"""
function coefficients_and_rhs_terms(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    ivz = cell_indices.ivz
    xstpavg = build_data.grid.xstpavg
    ystpavg = build_data.grid.ystpavg
    zstpavg = build_data.grid.zstpavg
    pscale = build_data.pscale
    bintern_zone = build_data.bc.bintern_zone
    bintern_velocity = build_data.bc.bintern_velocity
    Lii = build_data.system_vectors.Lii
    Ljj = build_data.system_vectors.Ljj
    Li = build_data.system_vectors.Li
    Lj = build_data.system_vectors.Lj
    Lv = build_data.system_vectors.Lv
    R = build_data.rhs.R

    kbond = 2.0 * pscale / (xstpavg + ystpavg + zstpavg)
    Lv[inz] = kbond
    Li[inz] = ivz
    Lj[inz] = ivz
    Lii[inz] = i
    Ljj[inz] = j
    inz = inz + 1
    if outside_prescribed_velocity_zone_z_stokes_3d(i, j, k, bintern_zone)
        R[ivz] = 0.0
    else
        R[ivz] = kbond * bintern_velocity[3]
    end
    return inz
end

end # module
