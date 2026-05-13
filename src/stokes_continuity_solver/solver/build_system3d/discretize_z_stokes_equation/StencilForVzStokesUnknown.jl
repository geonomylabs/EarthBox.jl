module StencilForVzStokesUnknown

include("boundary/BoundaryCells.jl")
include("internal/InternalCells.jl")

import EarthBox.BuildSysTools: SystemVectors
import ..Predicates3D: found_internal_cell_z_stokes_3d
import ..Predicates3D: outside_prescribed_velocity_zone_z_stokes_3d
import ..BuildStructs: StokesBuildData
import ..BuildStructs: CellIndices

""" Apply z-Stokes stencil for vz unknown from basic grid cell (i, j, k).

The z-Stokes stencil is applied to define coefficients and the right-hand
side term for the large matrix row ivz associated with the vz unknown for
the current grid cell with upper-left-front node (i, j, k). The vz unknown
is located in the center of the back z-perpendicular face of the basic
grid cell.

See the docstring of `StokesBuildManager3D` for the discretization of the
3D z-Stokes equation. This module is new in 3D and has no 2D analogue.

# Arguments
- `inz::Int`: Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).

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
    if found_internal_cell_z_stokes_3d(k, build_data.grid.znum) &&
       outside_prescribed_velocity_zone_z_stokes_3d(
           i, j, k, build_data.bc.bintern_zone)
        inz = InternalCells.coefficients_and_rhs_terms(inz, cell_indices, build_data)
    else
        inz = BoundaryCells.coefficients_and_rhs_terms(inz, cell_indices, build_data)
    end
    return inz
end

end # module
