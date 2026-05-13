module StencilForPContinuityUnknown

include("boundary/BoundaryCells.jl")
include("internal/InternalCells.jl")

import EarthBox.BuildSysTools: SystemVectors
import ..Predicates3D: found_internal_cell_continuity_3d
import ..BuildStructs: StokesBuildData
import ..BuildStructs: CellIndices

""" Apply continuity stencil for p unknown of basic grid cell (i, j, k).

See the docstring of `StokesBuildManager3D` for the discretization of the
3D continuity equation. The 3D continuity stencil samples six velocity
faces (vxL, vxR, vyU, vyD, vzF, vzB) of the basic cell with upper-left-
front node (i, j, k).

# Arguments
- `inz::Int`:
    - Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`:
    - Cell index information for the current cell.

- `build_data::StokesBuildData`:
    - Build data containing grid information, boundary conditions and
    arrays.

# Returns
- `inz::Int`: Updated index of non-zero matrix arrays.
"""
function apply_stencil(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    R = build_data.rhs.R
    RC = build_data.rhs.RC
    pscale = build_data.pscale
    ipr = cell_indices.ipr
    xnum = build_data.grid.xnum
    ynum = build_data.grid.ynum
    znum = build_data.grid.znum
    pressure_bc_mode = build_data.bc.pressure_bc_mode

    if found_internal_cell_continuity_3d(i, j, k, xnum, ynum, znum, pressure_bc_mode)
        R[ipr] = RC[i, j, k] * pscale
        inz = InternalCells.coefficients_and_rhs_terms(inz, cell_indices, build_data)
    else
        inz = BoundaryCells.coefficients_and_rhs_terms(inz, cell_indices, build_data)
    end
    return inz
end

end # module
