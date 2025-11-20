module InternalCells

include("stencil_nodes/CentralVelocityYNode.jl")
include("stencil_nodes/LeftVelocityYNode.jl")
include("stencil_nodes/RightVelocityYNode.jl")
include("stencil_nodes/UpperVelocityYNode.jl")
include("stencil_nodes/BottomVelocityYNode.jl")
include("stencil_nodes/UpperAndLowerPressureNodes.jl")
include("stencil_nodes/UpperLeftVelocityXNode.jl")
include("stencil_nodes/UpperRightVelocityXNode.jl")
include("stencil_nodes/BottomLeftVelocityXNode.jl")
include("stencil_nodes/BottomRightVelocityXNode.jl")

import EarthBox.BuildSysTools: SystemVectors
import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices

""" Calculate coefficients and rhs for vy unknowns in internal cells.

Internal cells are those with upper-left nodes located two steps
from right boundary or outside of prescribed velocity zones.

Coefficients and right-hand side term are calculated for the large matrix
row associated with the vy unknown for the current cell with upper-left
node (i,j). Coefficients are calculated for all unknown y-Stokes stencil
nodes and are applied to the appropriate matrix-row location associated
with vy. Coefficients are also updated for boundary conditions when
necessary. The right-hand-side term for the vy row is also calculated
taking boundary conditions into account.

The function applies only to basic cells with upper left node located at
least two steps from the bottom boundary of the basic grid. Consider the
following visual:

          o i=ynum-3
          |
          |
          |
          x i=ynum-2
          |
          |
          |
         LB i=ynum-1

where o refers to a node for which this function is applicable, x refers to
a non-applicable basic node and LB is a basic node located at lower
boundary of the basic grid.

See Figure 3 of the docstring in stokes_build.jl for details on the
y-Stokes stencil.

# Arguments
- `inz::Int`: Index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).

- `cell_indices::CellIndices`: 
    - Cell index information for the current cell.

- `build_data::StokesBuildData`: 
    - Build data.

# Returns
- `inz::Int`: Updated index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).
"""
function coefficients_and_rhs_terms(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    (
        drhody_gy_dt, drhodx_gy_dt
    ) = calculate_stabilization_terms(cell_indices, build_data)
    inz, inz_c = CentralVelocityYNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data, drhody_gy_dt)
    inz = UpperVelocityYNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = BottomVelocityYNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = LeftVelocityYNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = RightVelocityYNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz, inz_UL, inz_UR = UpperLeftVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data, drhodx_gy_dt)
    inz, inz_LL, inz_LR = BottomLeftVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data, drhodx_gy_dt)
    inz = UpperRightVelocityXNode.coefficients_and_rhs_term(
        inz, inz_UL, inz_UR, cell_indices, build_data, drhodx_gy_dt)
    inz = BottomRightVelocityXNode.coefficients_and_rhs_term(
        inz, inz_LL, inz_LR, cell_indices, build_data, drhodx_gy_dt)
    inz = UpperAndLowerPressureNodes.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    return inz
end

""" Calculate stabilization terms for vy unknowns in internal cells.

# Arguments
- `cell_indices::CellIndices`: Cell index information
- `build_data::StokesBuildData`: Build data

# Returns
- `Tuple{Float64,Float64}`: Tuple of (drhody_gy_dt, drhodx_gy_dt)
"""
function calculate_stabilization_terms(
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Tuple{Float64,Float64}
    if build_data.iuse_interface_stabilization == 1
        drhody = calculate_density_gradient_y(cell_indices, build_data)
        drhody_gy_dt = drhody * build_data.gravity_y * build_data.timestep

        drhodx = calculate_density_gradient_x(cell_indices, build_data)
        drhodx_gy_dt = drhodx * build_data.gravity_y * build_data.timestep
    else
        drhody_gy_dt = 0.0
        drhodx_gy_dt = 0.0
    end
    return drhody_gy_dt, drhodx_gy_dt
end

""" Calculate density gradient in y-direction.

# Arguments
- `cell_indices::CellIndices`: Cell index information
- `build_data::StokesBuildData`: Build data

# Returns
- `Float64`: Density gradient in y-direction
"""
function calculate_density_gradient_y(
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Float64
    i, j = cell_indices.i, cell_indices.j
    rho = build_data.rho1_vy
    ystp = build_data.grid.ystp

    drhody = (
        (rho[i+2, j+1] - rho[i+1, j+1])/ystp[i+1] +
        (rho[i+1, j+1] - rho[i, j+1])/ystp[i]
    )/2.0
    return drhody
end

""" Calculate density gradient in x-direction.

# Arguments
- `cell_indices::CellIndices`: Cell index information
- `build_data::StokesBuildData`: Build data

# Returns
- `Float64`: Density gradient in x-direction
"""
function calculate_density_gradient_x(
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Float64
    i = cell_indices.i
    j = cell_indices.j
    rho = build_data.rho1_vy
    xstpc = build_data.grid.xstpc

    drhodx = (
        (rho[i+1, j+2] - rho[i+1, j+1])/xstpc[j+1] +
        (rho[i+1, j+1] - rho[i+1, j])/xstpc[j]
    )/2.0
    return drhodx
end

end # module 