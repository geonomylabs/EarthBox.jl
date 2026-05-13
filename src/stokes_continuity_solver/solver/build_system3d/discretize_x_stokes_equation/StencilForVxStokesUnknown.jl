module StencilForVxStokesUnknown

include("boundary/BoundaryCells.jl")
include("internal/InternalCells.jl")

import EarthBox.BuildSysTools: SystemVectors
import ..Predicates3D: found_internal_cell_x_stokes_3d
import ..Predicates3D: outside_prescribed_velocity_zone_x_stokes_3d
import ..BuildStructs: StokesBuildData
import ..BuildStructs: CellIndices

""" Apply x-Stokes stencil for vx unknown from basic grid cell (i, j, k).

The x-Stokes stencil is applied to define coefficients and the right-hand
side term for the large matrix row ivx associated with the vx unknown for
the current grid cell with upper-left-front node (i, j, k). The vx unknown
is located in the center of the right x-perpendicular face of the basic
grid cell.

See the docstring of `StokesBuildManager3D` for the discretization of the
3D x-Stokes equation.

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
function apply_stencil(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    if found_internal_cell_x_stokes_3d(j, build_data.grid.xnum) &&
       outside_prescribed_velocity_zone_x_stokes_3d(
           i, j, k, build_data.bc.bintern_zone)
        inz = InternalCells.coefficients_and_rhs_terms(inz, cell_indices, build_data)
    else
        inz = BoundaryCells.coefficients_and_rhs_terms(inz, cell_indices, build_data)
    end
    return inz
end

end # module
