module InternalCells

include("stencil_nodes/LeftVelocityXNode.jl")
include("stencil_nodes/RightVelocityXNode.jl")
include("stencil_nodes/UpperVelocityYNode.jl")
include("stencil_nodes/BottomVelocityYNode.jl")
include("stencil_nodes/FrontVelocityZNode.jl")
include("stencil_nodes/BackVelocityZNode.jl")

import EarthBox.BuildSysTools: SystemVectors
import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices

""" Calculate coefficients and rhs for p unknowns in 3D internal cells.

Internal cells are basic cells whose pressure is not pinned by the
selected `pressure_bc_mode`. The 3D continuity stencil samples six
velocity faces of the cell and folds in boundary conditions for vx, vy
and vz when their face neighbours are at the model domain boundary.

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
    inz = LeftVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = RightVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = UpperVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = BottomVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = FrontVelocityZNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = BackVelocityZNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    return inz
end

end # module
