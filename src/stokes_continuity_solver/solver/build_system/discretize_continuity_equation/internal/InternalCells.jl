module InternalCells

include("stencil_nodes/LeftVelocityXNode.jl")
include("stencil_nodes/RightVelocityXNode.jl")
include("stencil_nodes/UpperVelocityYNode.jl")
include("stencil_nodes/BottomVelocityYNode.jl")

import EarthBox.BuildSysTools: SystemVectors
import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices

""" Calculate coefficients and rhs for p unknowns in internal cells.

Internal cells depend on pressure boundary condition type (pfirst_mode).
If pfirst_mode = 0, then the upper left corner cell is the boundary cell
while all others are internal cells for pressure. If pfirst_mode = 1 then
the top and bottom rows of cells are boundary cells while all others are
internal. If pfirst_mode = 2 the side columns of cells are boundary cells
while all others are internal.

The approach described in Figure 4 of the main module docstring is used to
discretize the continuity equation.

# Arguments
- `inz::Int`: Index of non-zero
             matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`: 
    - Cell index information for the current cell.

- `build_data::StokesBuildData`: 
    - Build data containing grid information, boundary conditions, and arrays.

# Returns
- `inz::Int`: Updated Index of non-zero
             matrix arrays (Lii, Ljj, Li, Lj and Lv).
"""
function coefficients_and_rhs_terms(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    inz = LeftVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = RightVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = UpperVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = BottomVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    return inz
end

end # module 