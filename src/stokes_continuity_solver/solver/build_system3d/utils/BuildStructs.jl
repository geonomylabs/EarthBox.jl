module BuildStructs

import Printf: @printf
import EarthBox.Arrays: ArrayUtils
import EarthBox.ModelDataContainer: ModelData
import EarthBox.BuildSysTools: SystemVectors

"""
CellIndices is a struct that contains the basic-grid indices and the
discretized large-matrix row indices for the four unknowns (vx, vy, vz, p')
associated with a 3D basic grid cell.

# Attributes
- `i::Int64`:
    - Basic grid node index for y-direction.

- `j::Int64`:
    - Basic grid node index for x-direction.

- `k::Int64`:
    - Basic grid node index for z-direction.

- `cell_index::Int64`:
    - Global basic grid cell index associated with (i, j, k).

- `ivx::Int64`:
    - Discretized large matrix index for the x-component of velocity vx
    associated with the current basic grid cell. The vx node is located in
    the center of the right-hand x-perpendicular face of the basic grid cell.

- `ivy::Int64`:
    - Discretized large matrix index for the y-component of velocity vy
    associated with basic grid cell index cell_index. The vy node is
    located in the center of the bottom y-perpendicular face of the basic
    grid cell.

- `ivz::Int64`:
    - Discretized large matrix index for the z-component of velocity vz
    associated with basic grid cell index cell_index. The vz node is
    located in the center of the front z-perpendicular face of the basic
    grid cell.

- `ipr::Int64`:
    - Discretized large matrix index for pressure pr. The pr node is located
    in the center of the current basic grid cell.
"""
struct CellIndices
    i::Int
    j::Int
    k::Int
    cell_index::Int
    ivx::Int
    ivy::Int
    ivz::Int
    ipr::Int

    """
    Calculate cell indices for a given i, j and k.

    # Arguments
    - `i::Int`: The i-index of the cell (y-direction).
    - `j::Int`: The j-index of the cell (x-direction).
    - `k::Int`: The k-index of the cell (z-direction).
    - `ynum::Int`: The number of basic grid nodes in the y-direction.
    - `xnum::Int`: The number of basic grid nodes in the x-direction.

    # Returns
    - `CellIndices`: A struct containing the cell indices.
    """
    function CellIndices(i::Int, j::Int, k::Int, ynum::Int, xnum::Int)::CellIndices
        cell_index = ((k - 1) * (xnum - 1) + (j - 1)) * (ynum - 1) + i
        ivx = cell_index * 4 - 3
        ivy = ivx + 1
        ivz = ivx + 2
        ipr = ivx + 3
        return new(i, j, k, cell_index, ivx, ivy, ivz, ipr)
    end
end

"""
GridData stores the basic-grid geometry, the spacing arrays needed by the
discretized stencils, and the matrix-row strides used to address neighbour
cells along the x and z axes.

# Attributes
- `xnum, ynum, znum::Int64`:
    - Number of basic grid nodes in each direction.

- `xstpavg, ystpavg, zstpavg::Float64`:
    - Average basic-grid spacings.

- `hshift_to_vxR::Int64`:
    - Matrix-row stride between a cell and its j+1 neighbour. Equal to
    (ynum-1)*4 in 3D.

- `dshift_to_vxF::Int64`:
    - Matrix-row stride between a cell and its k+1 neighbour. Equal to
    (xnum-1)*(ynum-1)*4 in 3D. New in 3D.

- `xstp, ystp, zstp::Vector{Float64}`:
    - Basic-grid spacings along each axis.

- `xstpc, ystpc, zstpc::Vector{Float64}`:
    - Spacings centered on basic nodes (used by staggered velocity faces).
"""
struct GridData
    xnum::Int64
    ynum::Int64
    znum::Int64
    xstpavg::Float64
    ystpavg::Float64
    zstpavg::Float64
    hshift_to_vxR::Int64
    dshift_to_vxF::Int64
    xstp::Vector{Float64}
    ystp::Vector{Float64}
    zstp::Vector{Float64}
    xstpc::Vector{Float64}
    ystpc::Vector{Float64}
    zstpc::Vector{Float64}
    function GridData(model::ModelData)
        xnum = model.grids.parameters.geometry.xnum.value
        ynum = model.grids.parameters.geometry.ynum.value
        znum = model.grids.parameters.geometry.znum.value
        hshift_to_vxR = (ynum - 1) * 4
        dshift_to_vxF = (xnum - 1) * (ynum - 1) * 4
        return new(
            xnum,
            ynum,
            znum,
            model.grids.parameters.geometry.xstpavg.value,
            model.grids.parameters.geometry.ystpavg.value,
            model.grids.parameters.geometry.zstpavg.value,
            hshift_to_vxR,
            dshift_to_vxF,
            model.grids.arrays.basic.xstp_b.array,
            model.grids.arrays.basic.ystp_b.array,
            model.grids.arrays.basic.zstp_b.array,
            model.grids.arrays.staggered_vy.xstp_vy.array,
            model.grids.arrays.staggered_vx.ystp_vx.array,
            model.grids.arrays.staggered_vx.zstp_vx.array
        )
    end
end

"""
BCData stores all boundary-condition arrays used by the 3D build system.

Each face array stores a (C, D) pair per face node, indexed as
`b<side><dir>[a, b, 1|2]` where the first two axes address the node on the
face and the last axis selects C (=1) or D (=2).

# Attributes
- `pressure_bc_mode::Int64`, `pressure_bc::Float64`:
    - Pressure boundary-condition mode and the prescribed pressure value.
    Mode 3 (front/back) is new in 3D.

- `btopx, btopy, btopz, bbottomx, bbottomy, bbottomz`:
    - Top/bottom face arrays (i = 1 and i = ynum boundaries).

- `bleftx, blefty, bleftz, brightx, brighty, brightz`:
    - Left/right face arrays (j = 1 and j = xnum boundaries).

- `bfrontx, bfronty, bfrontz, bbackx, bbacky, bbackz` (new in 3D):
    - Front/back face arrays (k = 1 and k = znum boundaries).

- `bintern_zone::Vector{Int64}`, `bintern_velocity::Vector{Float64}`:
    - Prescribed internal velocity zone bounds and the prescribed velocity
    components. See the module docstring of StokesBuildManager3D for the
    3D index conventions.
"""
struct BCData
    pressure_bc_mode::Int64
    pressure_bc::Float64
    btopx::Array{Float64,3}
    btopy::Array{Float64,3}
    btopz::Array{Float64,3}
    bbottomx::Array{Float64,3}
    bbottomy::Array{Float64,3}
    bbottomz::Array{Float64,3}
    bleftx::Array{Float64,3}
    blefty::Array{Float64,3}
    bleftz::Array{Float64,3}
    brightx::Array{Float64,3}
    brighty::Array{Float64,3}
    brightz::Array{Float64,3}
    bfrontx::Array{Float64,3}
    bfronty::Array{Float64,3}
    bfrontz::Array{Float64,3}
    bbackx::Array{Float64,3}
    bbacky::Array{Float64,3}
    bbackz::Array{Float64,3}
    bintern_zone::Vector{Int64}
    bintern_velocity::Vector{Float64}
    function BCData(model::ModelData)
        return new(
            model.bcs.parameters.bc_options.pressure_bc_mode.value,
            model.bcs.parameters.pressure.pressure_bc.value,
            model.bcs.arrays.vel_comp.btopx.array,
            model.bcs.arrays.vel_comp.btopy.array,
            model.bcs.arrays.vel_comp.btopz.array,
            model.bcs.arrays.vel_comp.bbottomx.array,
            model.bcs.arrays.vel_comp.bbottomy.array,
            model.bcs.arrays.vel_comp.bbottomz.array,
            model.bcs.arrays.vel_comp.bleftx.array,
            model.bcs.arrays.vel_comp.blefty.array,
            model.bcs.arrays.vel_comp.bleftz.array,
            model.bcs.arrays.vel_comp.brightx.array,
            model.bcs.arrays.vel_comp.brighty.array,
            model.bcs.arrays.vel_comp.brightz.array,
            model.bcs.arrays.vel_comp.bfrontx.array,
            model.bcs.arrays.vel_comp.bfronty.array,
            model.bcs.arrays.vel_comp.bfrontz.array,
            model.bcs.arrays.vel_comp.bbackx.array,
            model.bcs.arrays.vel_comp.bbacky.array,
            model.bcs.arrays.vel_comp.bbackz.array,
            model.bcs.arrays.internal.bintern_zone.array,
            model.bcs.arrays.internal.bintern_velocity.array
        )
    end
end

"""
RhsData stores the right-hand side terms of the discretized 3D system of
equations.

# Attributes
- `RC::Array{Float64,3}`, `RX::Array{Float64,3}`, `RY::Array{Float64,3}`,
  `RZ::Array{Float64,3}`:
    - Per-cell right-hand-side terms for the continuity, x-Stokes, y-Stokes
    and z-Stokes equations respectively.

- `R::Vector{Float64}`:
    - Flat right-hand side vector of the discretized system of equations.
"""
struct RhsData
    RC::Array{Float64,3}
    RX::Array{Float64,3}
    RY::Array{Float64,3}
    RZ::Array{Float64,3}
    R::Vector{Float64}
    function RhsData(model::ModelData)
        return new(
            model.stokes_continuity.arrays.rhs.RC1.array,
            model.stokes_continuity.arrays.rhs.RX1.array,
            model.stokes_continuity.arrays.rhs.RY1.array,
            model.stokes_continuity.arrays.rhs.RZ1.array,
            model.stokes_continuity.arrays.rhs.RHS.array,
        )
    end
end

"""
StokesBuildData groups all data needed by the 3D Stokes-continuity build
system. It is constructed once per build pass and threaded through every
stencil module.

The pressure scale `pscale` is recomputed locally with the 3D formula
`pscale = 3 * etan[1,1,1] / (xstpavg + ystpavg + zstpavg)` rather than read
from the upstream `ModelData` parameter, which still carries the 2D scaling.

Three shear viscosities are required, one per shear-stress component:
`etas_xy` on z-parallel basic edges, `etas_xz` on y-parallel basic edges and
`etas_yz` on x-parallel basic edges. The 2D `etas` field corresponds to
`etas_xy` in this 3D notation.

Density samples are kept per velocity face (`rho1_vx`, `rho1_vy`,
`rho1_vz`). Only `rho1_vy` participates in the gravity-driven right-hand
side and interface-stabilization terms today, but the other two are
included for the new x- and z-momentum equations.
"""
struct StokesBuildData
    timestep::Float64
    pscale::Float64
    gravity_x::Float64
    gravity_y::Float64
    gravity_z::Float64
    iuse_interface_stabilization::Int64
    grid::GridData
    bc::BCData
    rhs::RhsData
    system_vectors::SystemVectors
    etan::Array{Float64,3}
    etas_xy::Array{Float64,3}
    etas_xz::Array{Float64,3}
    etas_yz::Array{Float64,3}
    rho1_vx::Array{Float64,3}
    rho1_vy::Array{Float64,3}
    rho1_vz::Array{Float64,3}
    function StokesBuildData(model::ModelData)
        grid = GridData(model)
        etan = model.stokes_continuity.arrays.viscosity.etan0.array
        # 3D pressure scaling. The first interior pressure node samples
        # etan, which is consistent with the 2D convention applied to
        # etan[1, 1].
        pscale = 3.0 * etan[1, 1, 1] / (grid.xstpavg + grid.ystpavg + grid.zstpavg)
        return new(
            model.timestep.parameters.main_time_loop.timestep.value,
            pscale,
            model.gravity.parameters.gravity_x.value,
            model.gravity.parameters.gravity_y.value,
            model.gravity.parameters.gravity_z.value,
            model.stokes_continuity.parameters.build.iuse_interface_stabilization.value,
            grid,
            BCData(model),
            RhsData(model),
            model.stokes_continuity.parameters.build.system_vectors,
            etan,
            model.stokes_continuity.arrays.viscosity.etas_xy0.array,
            model.stokes_continuity.arrays.viscosity.etas_xz0.array,
            model.stokes_continuity.arrays.viscosity.etas_yz0.array,
            model.stokes_continuity.arrays.density.rho1_vx.array,
            model.stokes_continuity.arrays.density.rho1_vy.array,
            model.stokes_continuity.arrays.density.rho1_vz.array
        )
    end
end

end
