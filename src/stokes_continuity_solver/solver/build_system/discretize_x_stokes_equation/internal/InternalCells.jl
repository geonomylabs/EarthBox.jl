module InternalCells

import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices

include("stencil_nodes/CentralVelocityXNode.jl")
include("stencil_nodes/LeftVelocityXNode.jl")
include("stencil_nodes/RightVelocityXNode.jl")
include("stencil_nodes/UpperVelocityXNode.jl")
include("stencil_nodes/BottomVelocityXNode.jl")
include("stencil_nodes/UpperLeftVelocityYNode.jl")
include("stencil_nodes/UpperRightVelocityYNode.jl")
include("stencil_nodes/BottomLeftVelocityYNode.jl")
include("stencil_nodes/BottomRightVelocityYNode.jl")
include("stencil_nodes/LeftAndRightPressureNodes.jl")

""" Calculate coefficients and rhs for vx unknowns in internal cells.

Internal cells are those cells with upper-left nodes located two steps
from right boundary or outside of prescribed velocity zones.

Coefficients and right-hand side terms are calculated for the large matrix
row associated with the vx unknown for the current cell with upper-left
node (i,j). Coefficients are calculated for all unknown x-Stokes stencil
nodes and are applied to the appropriate matrix-row location associated
with vx. Coefficients are also updated for boundary conditions when
necessary. The right-hand-side term for the vx row is also calculated
taking boundary conditions into account.

See Figure 2 of the docstring in stokes_build.jl for details on
the x-Stokes stencil.

# Arguments
- `inz::Int`: Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`: 
    - Cell index information for the current cell.

- `build_data::StokesBuildData`: 
    - Build data.

# Returns
- `inz`: Updated Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).

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
    inz, inz_UL, inz_LL = UpperLeftVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz, inz_UR, inz_LR = UpperRightVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    inz = BottomLeftVelocityYNode.coefficients_and_rhs_term(
        inz, inz_UL, inz_LL, cell_indices, build_data)
    inz = BottomRightVelocityYNode.coefficients_and_rhs_term(
        inz, inz_UR, inz_LR, cell_indices, build_data)
    inz = LeftAndRightPressureNodes.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    return inz
end

end # module