module InternalCells

include("stencil_nodes/CentralVelocityZNode.jl")
include("stencil_nodes/FrontVelocityZNode.jl")
include("stencil_nodes/BackVelocityZNode.jl")
include("stencil_nodes/UpperVelocityZNode.jl")
include("stencil_nodes/BottomVelocityZNode.jl")
include("stencil_nodes/LeftVelocityZNode.jl")
include("stencil_nodes/RightVelocityZNode.jl")
include("stencil_nodes/FrontUpperVelocityYNode.jl")
include("stencil_nodes/FrontBottomVelocityYNode.jl")
include("stencil_nodes/BackUpperVelocityYNode.jl")
include("stencil_nodes/BackBottomVelocityYNode.jl")
include("stencil_nodes/FrontLeftVelocityXNode.jl")
include("stencil_nodes/FrontRightVelocityXNode.jl")
include("stencil_nodes/BackLeftVelocityXNode.jl")
include("stencil_nodes/BackRightVelocityXNode.jl")
include("stencil_nodes/FrontAndBackPressureNodes.jl")

import EarthBox.BuildSysTools: SystemVectors
import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices

""" Calculate coefficients and rhs for vz unknowns in 3D internal cells.

Internal cells are those with upper-left-front nodes located two basic
slabs from the back z-boundary or outside any prescribed-vz zone.

Coefficients and right-hand-side terms are calculated for the large matrix
row associated with the vz unknown for the current cell with upper-left-
front node (i, j, k). Coefficients are calculated for all unknown nodes
in the 3D z-Stokes stencil (six vz neighbours, four vy edge neighbours
for the yz-shear, four vx edge neighbours for the xz-shear, and the
front/back pressure pair). Coefficients are also updated for boundary
conditions when necessary. The right-hand side term for the vz row is
also calculated taking boundary conditions into account.

Gravity acts along y, so this stencil does not carry an interface-
stabilization correction; that correction lives only in the y-Stokes
stencil.

See the module docstring of `StokesBuildManager3D` for the 3D z-Stokes
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
    inz, inz_c = CentralVelocityZNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = FrontVelocityZNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = BackVelocityZNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = UpperVelocityZNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = BottomVelocityZNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = LeftVelocityZNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = RightVelocityZNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz, inz_FU, inz_FD = FrontUpperVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz, inz_BU, inz_BD = BackUpperVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = FrontBottomVelocityYNode.coefficients_and_rhs_term(
        inz, inz_FU, inz_FD, cell_indices, build_data)
    inz = BackBottomVelocityYNode.coefficients_and_rhs_term(
        inz, inz_BU, inz_BD, cell_indices, build_data)
    inz, inz_FL, inz_FR = FrontLeftVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz, inz_BL, inz_BR = BackLeftVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = FrontRightVelocityXNode.coefficients_and_rhs_term(
        inz, inz_FL, inz_FR, cell_indices, build_data)
    inz = BackRightVelocityXNode.coefficients_and_rhs_term(
        inz, inz_BL, inz_BR, cell_indices, build_data)
    inz = FrontAndBackPressureNodes.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    return inz
end

end # module
