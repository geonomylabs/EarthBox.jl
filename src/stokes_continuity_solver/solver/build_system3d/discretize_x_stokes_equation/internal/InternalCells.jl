module InternalCells

import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices

include("stencil_nodes/CentralVelocityXNode.jl")
include("stencil_nodes/LeftVelocityXNode.jl")
include("stencil_nodes/RightVelocityXNode.jl")
include("stencil_nodes/UpperVelocityXNode.jl")
include("stencil_nodes/BottomVelocityXNode.jl")
include("stencil_nodes/FrontVelocityXNode.jl")
include("stencil_nodes/BackVelocityXNode.jl")
include("stencil_nodes/UpperLeftVelocityYNode.jl")
include("stencil_nodes/UpperRightVelocityYNode.jl")
include("stencil_nodes/BottomLeftVelocityYNode.jl")
include("stencil_nodes/BottomRightVelocityYNode.jl")
include("stencil_nodes/FrontLeftVelocityZNode.jl")
include("stencil_nodes/FrontRightVelocityZNode.jl")
include("stencil_nodes/BackLeftVelocityZNode.jl")
include("stencil_nodes/BackRightVelocityZNode.jl")
include("stencil_nodes/LeftAndRightPressureNodes.jl")

""" Calculate coefficients and rhs for vx unknowns in 3D internal cells.

Internal cells are those with upper-left-front nodes located two basic
columns from the right x-boundary or outside any prescribed-vx zone.

Coefficients and right-hand-side terms are calculated for the large matrix
row associated with the vx unknown for the current cell with upper-left-
front node (i, j, k). Coefficients are calculated for all unknown nodes in
the 3D x-Stokes stencil (six vx neighbours, four vy edge neighbours, four
vz edge neighbours and the left/right pressure pair). Coefficients are
also updated for boundary conditions when necessary. The right-hand-side
term for the vx row is also calculated taking boundary conditions into
account.

See the module docstring of `StokesBuildManager3D` for the 3D x-Stokes
stencil.

# Arguments
- `inz::Int`: Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).

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
    inz, inz_c = CentralVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = LeftVelocityXNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = RightVelocityXNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = UpperVelocityXNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = BottomVelocityXNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = FrontVelocityXNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = BackVelocityXNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz, inz_UL, inz_LL = UpperLeftVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz, inz_UR, inz_LR = UpperRightVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = BottomLeftVelocityYNode.coefficients_and_rhs_term(
        inz, inz_UL, inz_LL, cell_indices, build_data)
    inz = BottomRightVelocityYNode.coefficients_and_rhs_term(
        inz, inz_UR, inz_LR, cell_indices, build_data)
    inz, inz_FL, inz_BL = FrontLeftVelocityZNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz, inz_FR, inz_BR = FrontRightVelocityZNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = BackLeftVelocityZNode.coefficients_and_rhs_term(
        inz, inz_FL, inz_BL, cell_indices, build_data)
    inz = BackRightVelocityZNode.coefficients_and_rhs_term(
        inz, inz_FR, inz_BR, cell_indices, build_data)
    inz = LeftAndRightPressureNodes.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    return inz
end

end # module
