module StencilForVyStokesUnknown

include("boundary/BoundaryCells.jl")
include("internal/InternalCells.jl")

import EarthBox.BuildSysTools: SystemVectors
import EarthBox.Domain: found_internal_cell_y_stokes
import EarthBox.Domain: outside_prescribed_velocity_zone_y_stokes
import ..BuildStructs: StokesBuildData
import ..BuildStructs: CellIndices

""" Apply y-Stokes stencil for vy unknown from basic grid cell (i,j).

The y-Stokes stencil is applied to defined coefficients and right-hand side
term for the large matrix row ivy associated with the vy unknown for the
current grid cell with upper left node (i,j). The vy unknown is located
along the middle of the bottom boundary of the basic grid cell.

The approach described in Figure 3 of the docstring in stokes_build.jl is
used to discretize the y-Stokes equation.

# Arguments
- `inz::Int`: 
    - Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`: 
    - Cell index information for the current cell.

- `build_data::StokesBuildData`: 
    - Build data.

# Returns
- `inz::Int`: Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).
"""
function apply_stencil(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    if found_internal_cell_y_stokes(i, build_data.grid.ynum) &&
       outside_prescribed_velocity_zone_y_stokes(i, j, build_data.bc.bintern_zone)
        inz = InternalCells.coefficients_and_rhs_terms(inz, cell_indices, build_data)
    else
        inz = BoundaryCells.coefficients_and_rhs_terms(inz, cell_indices, build_data)
    end
    return inz
end

end # module 