module BuildStructs

import Printf: @printf
import Plots
import EarthBox.Arrays: ArrayUtils
import EarthBox.ModelDataContainer: ModelData
import EarthBox.BuildSysTools: SystemVectors
import ...GlobalIndices

"""
CellIndices is a struct that contains the.

# Attributes
- `i::Int64`: 
    - Basic grid node index for y-direction.

- `j::Int64`: 
    - Basic grid node index for x-direction.

- `cell_index::Int64`: 
    - Global basic grid cell index associated with (i, j) (see numbers in upper 
    right-hand corner of grid cells in Figure 1 in StokesBuildManager.jl).

- `ivx::Int64`: 
    - Discretized large matrix index for the x-component of velocity vx
    associated with the current basic grid cell. The vx node is located in
    the middle of the right-hand boundary of the basic grid cell.

- `ivy::Int64`: 
    - Discretized large matrix index for the y-component of velocity vy
    associated with basic grid cell index cell_index. The vy node is
    located in the middle of the bottom boundary of the basic grid cell.

- `ipr::Int64`: 
    - Discretized large matrix index for pressure pr. The pr node is located
    in the middle of the current basic grid cell.
"""
struct CellIndices
    i::Int
    j::Int
    cell_index::Int
    ivx::Int
    ivy::Int
    ipr::Int

    """
    Calculate cell indices for a given i and j.

    # Arguments
    - `i::Int`: The i-index of the cell (y-direction).
    - `j::Int`: The j-index of the cell (x-direction).
    - `ynum::Int`: The number basic grid nodes in the y-direction.

    # Returns
    - `CellIndices`: A struct containing the cell indices.
    """
    function CellIndices(i::Int, j::Int, ynum::Int)::CellIndices
        cell_index = GlobalIndices.get_global_basic_cell_index(i, j, ynum)
        ivx = GlobalIndices.get_global_ivx_unknown_index(cell_index)
        ivy = GlobalIndices.get_global_ivy_unknown_index(ivx)
        ipr = GlobalIndices.get_global_ipr_unknown_index(ivx)
        return new(i, j, cell_index, ivx, ivy, ipr)
    end
end

struct GridData
    xnum::Int64
    ynum::Int64
    xstpavg::Float64
    ystpavg::Float64
    hshift_to_vxR::Float64
    xstp::Vector{Float64}
    ystp::Vector{Float64}
    xstpc::Vector{Float64}
    ystpc::Vector{Float64}
    function GridData(model::ModelData)
        return new(
            model.grids.parameters.geometry.xnum.value,
            model.grids.parameters.geometry.ynum.value,
            model.grids.parameters.geometry.xstpavg.value,
            model.grids.parameters.geometry.ystpavg.value,
            model.stokes_continuity.parameters.build.hshift_to_vxR.value,
            model.grids.arrays.basic.xstp_b.array,
            model.grids.arrays.basic.ystp_b.array,
            model.grids.arrays.staggered_vy.xstp_vy.array,
            model.grids.arrays.staggered_vx.ystp_vx.array
        )
    end
end

struct BCData
    pressure_bc_mode::Int64
    pressure_bc::Float64
    btopx::Matrix{Float64}
    btopy::Matrix{Float64}
    bbottomx::Matrix{Float64}
    bbottomy::Matrix{Float64}
    bleftx::Matrix{Float64}
    blefty::Matrix{Float64}
    brightx::Matrix{Float64}
    brighty::Matrix{Float64}
    bintern_zone::Vector{Int64}
    bintern_velocity::Vector{Float64}
    function BCData(model::ModelData)
        return new(
            model.bcs.parameters.bc_options.pressure_bc_mode.value,
            model.bcs.parameters.pressure.pressure_bc.value,
            model.bcs.arrays.vel_comp.btopx.array,
            model.bcs.arrays.vel_comp.btopy.array,
            model.bcs.arrays.vel_comp.bbottomx.array,
            model.bcs.arrays.vel_comp.bbottomy.array,
            model.bcs.arrays.vel_comp.bleftx.array,
            model.bcs.arrays.vel_comp.blefty.array,
            model.bcs.arrays.vel_comp.brightx.array,
            model.bcs.arrays.vel_comp.brighty.array,
            model.bcs.arrays.internal.bintern_zone.array,
            model.bcs.arrays.internal.bintern_velocity.array
        )
    end
end

"""
RhsData is a struct that contains the right-hand side of the discretized system 
of equations.

# Attributes
- `RC::Matrix{Float64}`: 

- `RX::Matrix{Float64}`: 

- `RY::Matrix{Float64}`: 

- `R::Vector{Float64}`: 
    - Right-hand side array of the discretized system of equations.
"""
struct RhsData
    RC::Matrix{Float64}
    RX::Matrix{Float64}
    RY::Matrix{Float64}
    R::Vector{Float64}
    function RhsData(model::ModelData)
        return new(
            model.stokes_continuity.arrays.rhs.RC1.array,
            model.stokes_continuity.arrays.rhs.RX1.array,
            model.stokes_continuity.arrays.rhs.RY1.array,
            model.stokes_continuity.arrays.rhs.RHS.array,
        )
    end
end

struct StokesBuildData
    timestep::Float64
    pscale::Float64
    gravity_y::Float64
    iuse_interface_stabilization::Int64
    grid::GridData
    bc::BCData
    rhs::RhsData
    system_vectors::SystemVectors
    etan::Matrix{Float64}
    etas::Matrix{Float64}
    rho1_vy::Matrix{Float64}
    function StokesBuildData(model::ModelData)
        return new(
            model.timestep.parameters.main_time_loop.timestep.value,
            model.stokes_continuity.parameters.build.pscale.value,
            model.gravity.parameters.gravity_y.value,
            model.stokes_continuity.parameters.build.iuse_interface_stabilization.value,
            GridData(model),
            BCData(model),
            RhsData(model),
            SystemVectors(model.stokes_continuity.parameters.build.nonzero_max.value),
            model.stokes_continuity.arrays.viscosity.etan0.array,
            model.stokes_continuity.arrays.viscosity.etas0.array,
            model.stokes_continuity.arrays.density.rho1_vy.array
        )
    end
end

end

