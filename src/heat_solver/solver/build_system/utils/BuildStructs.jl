module BuildStructs

import EarthBox.ModelDataContainer: ModelData
import EarthBox.BuildSysTools: SystemVectors

struct CellIndices
    i::Int
    j::Int
    itk::Int
end
    
struct GridData
    xnum::Int64
    ynum::Int64
    xstp::Vector{Float64}
    ystp::Vector{Float64}
    xstpc::Vector{Float64}
    ystpc::Vector{Float64}
    function GridData(model::ModelData)
        return new(
            model.grids.parameters.geometry.xnum.value,
            model.grids.parameters.geometry.ynum.value,
            model.grids.arrays.basic.xstp_b.array,
            model.grids.arrays.basic.ystp_b.array,
            model.grids.arrays.pressure.xstp_pr.array,
            model.grids.arrays.pressure.ystp_pr.array
        )
    end
end

struct BCData
    btop::Matrix{Float64}
    bbottom::Matrix{Float64}
    bleft::Matrix{Float64}
    bright::Matrix{Float64}
    function BCData(model::ModelData)
        return new(
            model.bcs.arrays.temperature.btopt.array,
            model.bcs.arrays.temperature.bbottomt.array,
            model.bcs.arrays.temperature.bleftt.array,
            model.bcs.arrays.temperature.brightt.array
        )
    end
end

struct RhsData
    RT::Matrix{Float64}
    R::Vector{Float64}
    function RhsData(model::ModelData)
        return new(
            model.heat_equation.arrays.rhs.RT1.array,
            model.heat_equation.arrays.rhs.RHSheat.array
        )
    end
end

struct HeatBuildData
    # Time step in seconds used to solve the transient heat conduction equation.
    # The current heat loop time step timestep_heat is used. The heat loop time step
    # is initialized to the model time step at he beginning of the heat loop.
    timestep::Float64
    grid::GridData
    bc::BCData
    rhs::RhsData
    system_vectors::SystemVectors
    tk::Matrix{Float64}
    rhocp::Matrix{Float64}
    kt::Matrix{Float64}
    function HeatBuildData(model::ModelData)
        return new(
            model.timestep.parameters.thermal_loop.timestep_heat.value,
            GridData(model),
            BCData(model),
            RhsData(model),
            SystemVectors(model.heat_equation.parameters.build.nonzero_max_heat.value),
            model.heat_equation.arrays.temperature.tk0.array,
            model.heat_equation.arrays.rhocp.rhocp1.array,
            model.heat_equation.arrays.thermal_conductivity.kt1.array
        )
    end
end

"""
Define scalar information at stencil nodes where C refers to central, U refers 
to upper, D refers to down, L refers to left and R refers to right. U, D, L, R 
are defined relative to the central node.
"""
struct HeatStencilData
    tkC::Float64
    rhocpC::Float64
    ktU::Float64
    ktC::Float64
    ktD::Float64
    ktL::Float64
    ktR::Float64
    dxL::Float64
    dxR::Float64
    dxC::Float64
    dyU::Float64
    dyD::Float64
    dyC::Float64
    timestep::Float64
    ynum::Int
    function HeatStencilData(build_data::HeatBuildData, cell_indices::CellIndices)
        i = cell_indices.i
        j = cell_indices.j
        return new(
            build_data.tk[i, j],
            build_data.rhocp[i, j],
            build_data.kt[i-1, j],
            build_data.kt[i, j],
            build_data.kt[i+1, j],
            build_data.kt[i, j-1],
            build_data.kt[i, j+1],
            build_data.grid.xstp[j-1],
            build_data.grid.xstp[j],
            build_data.grid.xstpc[j-1],
            build_data.grid.ystp[i-1],
            build_data.grid.ystp[i],
            build_data.grid.ystpc[i-1],
            build_data.timestep,
            build_data.grid.ynum
        )
    end
end

end # module