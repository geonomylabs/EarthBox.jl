module InternalCells

include("stencil_nodes/CentralVelocityYNode.jl")
include("stencil_nodes/LeftVelocityYNode.jl")
include("stencil_nodes/RightVelocityYNode.jl")
include("stencil_nodes/UpperVelocityYNode.jl")
include("stencil_nodes/BottomVelocityYNode.jl")
include("stencil_nodes/FrontVelocityYNode.jl")
include("stencil_nodes/BackVelocityYNode.jl")
include("stencil_nodes/UpperAndLowerPressureNodes.jl")
include("stencil_nodes/UpperLeftVelocityXNode.jl")
include("stencil_nodes/UpperRightVelocityXNode.jl")
include("stencil_nodes/BottomLeftVelocityXNode.jl")
include("stencil_nodes/BottomRightVelocityXNode.jl")
include("stencil_nodes/FrontUpperVelocityZNode.jl")
include("stencil_nodes/FrontBottomVelocityZNode.jl")
include("stencil_nodes/BackUpperVelocityZNode.jl")
include("stencil_nodes/BackBottomVelocityZNode.jl")

import EarthBox.BuildSysTools: SystemVectors
import ...BuildStructs: StokesBuildData
import ...BuildStructs: CellIndices

""" Calculate coefficients and rhs for vy unknowns in 3D internal cells.

Internal cells are those with upper-left-front nodes located two basic
rows from the lower y-boundary or outside any prescribed-vy zone.

Coefficients and right-hand-side terms are calculated for the large matrix
row associated with the vy unknown for the current cell with upper-left-
front node (i, j, k). Coefficients are calculated for all unknown nodes
in the 3D y-Stokes stencil (six vy neighbours, four vx edge neighbours,
four vz edge neighbours and the upper/lower pressure pair). Coefficients
are also updated for boundary conditions when necessary. The right-hand
side term for the vy row is also calculated taking boundary conditions
into account.

Gravity acts along y, so the interface-stabilization correction
(`drhody_gy_dt` on the diagonal, `drhodx_gy_dt` on the four vx edge
neighbours, `drhodz_gy_dt` on the four vz edge neighbours) is applied
here.

See the module docstring of `StokesBuildManager3D` for the 3D y-Stokes
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
    drhody_gy_dt, drhodx_gy_dt, drhodz_gy_dt = calculate_stabilization_terms(
        cell_indices, build_data)
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
    inz = FrontVelocityYNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz = BackVelocityYNode.coefficients_and_rhs_term(
        inz, inz_c, cell_indices, build_data)
    inz, inz_UL, inz_UR = UpperLeftVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data, drhodx_gy_dt)
    inz, inz_LL, inz_LR = BottomLeftVelocityXNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data, drhodx_gy_dt)
    inz = UpperRightVelocityXNode.coefficients_and_rhs_term(
        inz, inz_UL, inz_UR, cell_indices, build_data, drhodx_gy_dt)
    inz = BottomRightVelocityXNode.coefficients_and_rhs_term(
        inz, inz_LL, inz_LR, cell_indices, build_data, drhodx_gy_dt)
    inz, inz_FU, inz_BU = FrontUpperVelocityZNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data, drhodz_gy_dt)
    inz, inz_FD, inz_BD = FrontBottomVelocityZNode.coefficients_and_rhs_term(
        inz, cell_indices, build_data, drhodz_gy_dt)
    inz = BackUpperVelocityZNode.coefficients_and_rhs_term(
        inz, inz_FU, inz_BU, cell_indices, build_data, drhodz_gy_dt)
    inz = BackBottomVelocityZNode.coefficients_and_rhs_term(
        inz, inz_FD, inz_BD, cell_indices, build_data, drhodz_gy_dt)
    inz = UpperAndLowerPressureNodes.coefficients_and_rhs_term(
        inz, cell_indices, build_data)
    return inz
end

""" Calculate interface-stabilization terms for vy unknowns in internal cells.

Three centered density gradients are computed (drhody, drhodx, drhodz)
and scaled by gravity_y * timestep. The y-direction gradient feeds the
diagonal vy coefficient, the x-direction gradient feeds the four xy-shear
vx edge coefficients and the z-direction gradient feeds the four yz-shear
vz edge coefficients.

# Returns
- `Tuple{Float64,Float64,Float64}`: (drhody_gy_dt, drhodx_gy_dt, drhodz_gy_dt)
"""
function calculate_stabilization_terms(
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Tuple{Float64,Float64,Float64}
    if build_data.iuse_interface_stabilization == 1
        drhody = calculate_density_gradient_y(cell_indices, build_data)
        drhody_gy_dt = drhody * build_data.gravity_y * build_data.timestep

        drhodx = calculate_density_gradient_x(cell_indices, build_data)
        drhodx_gy_dt = drhodx * build_data.gravity_y * build_data.timestep

        drhodz = calculate_density_gradient_z(cell_indices, build_data)
        drhodz_gy_dt = drhodz * build_data.gravity_y * build_data.timestep
    else
        drhody_gy_dt = 0.0
        drhodx_gy_dt = 0.0
        drhodz_gy_dt = 0.0
    end
    return drhody_gy_dt, drhodx_gy_dt, drhodz_gy_dt
end

""" Calculate density gradient in y-direction at vyC. """
function calculate_density_gradient_y(
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Float64
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    rho = build_data.rho1_vy
    ystp = build_data.grid.ystp

    drhody = (
        (rho[i+2, j+1, k+1] - rho[i+1, j+1, k+1]) / ystp[i+1] +
        (rho[i+1, j+1, k+1] - rho[i, j+1, k+1]) / ystp[i]
    ) / 2.0
    return drhody
end

""" Calculate density gradient in x-direction at vyC. """
function calculate_density_gradient_x(
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Float64
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    rho = build_data.rho1_vy
    xstpc = build_data.grid.xstpc

    drhodx = (
        (rho[i+1, j+2, k+1] - rho[i+1, j+1, k+1]) / xstpc[j+1] +
        (rho[i+1, j+1, k+1] - rho[i+1, j, k+1]) / xstpc[j]
    ) / 2.0
    return drhodx
end

""" Calculate density gradient in z-direction at vyC. """
function calculate_density_gradient_z(
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Float64
    i = cell_indices.i
    j = cell_indices.j
    k = cell_indices.k
    rho = build_data.rho1_vy
    zstpc = build_data.grid.zstpc

    drhodz = (
        (rho[i+1, j+1, k+2] - rho[i+1, j+1, k+1]) / zstpc[k+1] +
        (rho[i+1, j+1, k+1] - rho[i+1, j+1, k]) / zstpc[k]
    ) / 2.0
    return drhodz
end

end # module
