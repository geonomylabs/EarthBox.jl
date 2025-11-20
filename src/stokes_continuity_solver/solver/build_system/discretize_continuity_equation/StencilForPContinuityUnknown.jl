module StencilForPContinuityUnknown

include("boundary/BoundaryCells.jl")
include("internal/InternalCells.jl")

import EarthBox.BuildSysTools: SystemVectors
import EarthBox.Domain: found_internal_cell_continuity
import ..BuildStructs: StokesBuildData
import ..BuildStructs: CellIndices

""" Apply continuity stencil for p unknown of basic grid cell (i,j).

The approach described in Figure 4 of the main module docstring is used to
discretize the continuity equation.

# Arguments
- `inz::Int`: 
    - Index of non-zero
      matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`: 
    - Cell index information for the current cell.

- `build_data::StokesBuildData`: 
    - Build data containing grid information, boundary conditions, and arrays.

# Returns
- `inz::Int`: Index of non-zero
             matrix arrays (Lii, Ljj, Li, Lj and Lv).
"""
function apply_stencil(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    R = build_data.rhs.R
    RC = build_data.rhs.RC
    pscale = build_data.pscale
    ipr = cell_indices.ipr
    xnum = build_data.grid.xnum
    ynum = build_data.grid.ynum
    pressure_bc_mode = build_data.bc.pressure_bc_mode

    if found_internal_cell_continuity(i, j, xnum, ynum, pressure_bc_mode)
        R[ipr] = RC[i, j] * pscale
        inz = InternalCells.coefficients_and_rhs_terms(inz, cell_indices, build_data)
    else
        inz = BoundaryCells.coefficients_and_rhs_terms(inz, cell_indices, build_data)
    end
    return inz
end

end # module 