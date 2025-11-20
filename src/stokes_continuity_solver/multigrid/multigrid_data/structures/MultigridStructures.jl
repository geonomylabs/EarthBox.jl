module MultigridStructures

import EarthBox.Arrays.ArrayTypes.ScalarArray2D: grid_array2D
import EarthBox.Arrays.ArrayTypes.ScalarArray3D: grid_array3D
import ...LevelManager
import ...LevelManager: LevelData, LevelData2d

Base.@kwdef mutable struct ResidualScalingFactors
    stokesscale::Float64 = 1.0
    continscale::Float64 = 1.0
end

Base.@kwdef mutable struct Counters
    # Counter for v-cycles.
    ivcycle::Int64=1
    # Counter for global updates for multi-multigrid.
    iglobal_update::Int64=1
    # V-cycle counter for viscosity jumps.
    iviscosity_jump::Int64=1
    # Counter of maximum viscosity oversteps for multi-multigrid. When the
    # number of oversteps is greater than 1 a rest occurs and the global
    # solution is updated.
    iviscosity_max_overstep::Int64=0
end

Base.@kwdef mutable struct ViscosityScalingParameters
    # Current minimum computational viscosity Pa.s
    viscosity_min_cur::Float64=1e20
    # Current maximum computational viscosity Pa.s
    viscosity_max_cur::Float64=1e20
    # Starting maximum computational viscosity Pa.s
    viscosity_max_start::Float64=1e20
    # Step size factor for minimum computational viscosity
    viscosity_min_factor::Float64=1.0
    # Step size factor for maximum computational viscosity
    viscosity_max_factor::Float64=3.3333
    # Original (starting) minimum viscosity of model Pa.s
    viscosity_min_o::Float64=1e20
    # Original (starting) maximum viscosity of model Pa.s
    viscosity_max_o::Float64=1e23
    # Number of viscosity jumps
    nviscosity_jumps::Float64=12.0
    # Number of v-cycles between viscosity jumps
    nvcycles_viscosity_jump::Int64=15
    # Maximum residual for 3D viscosity jumps
    resmax::Float64=1e-6
end

function ViscosityScalingParameters(vcycle::NamedTuple)::ViscosityScalingParameters
    return ViscosityScalingParameters(
        nviscosity_jumps=vcycle.nviscosity_jumps,
        nvcycles_viscosity_jump=vcycle.nvcycles_viscosity_jump
    )
end

mutable struct VcycleParameters
    # Boolean flag used to determine if multi-multigrid is used.
    use_multimulti::Bool
    # Number of v-cycles.
    nvcycles::Int64
    # Convergence criterion after viscosity jumps are completed.
    convergence_criterion::Float64
    # Boolean flag used to determine if plots are made during v-cycle calculations
    make_plots::Bool
end

function VcycleParameters(vcycle::NamedTuple)::VcycleParameters
    return VcycleParameters(
        vcycle.use_multimulti,
        vcycle.nvcycles,
        vcycle.convergence_criterion,
        vcycle.make_plots
    )
end

struct RelaxationParameters
    # Relaxation parameter for Stokes equation (old name = krelaxs)
    relax_stokes::Float64
    # Relaxation parameter for velocity equation (old name = vkoef)
    relax_velocity::Float64
    # Relaxation parameter for continuity equation (old name = krelaxc)
    relax_continuity::Float64
    # Relaxation parameter for pressure equation (old name = pkoef)
    relax_pressure::Float64
end

function RelaxationParameters(relaxation::NamedTuple)::RelaxationParameters
    return RelaxationParameters(
        relaxation.relax_stokes, relaxation.relax_velocity, 
        relaxation.relax_continuity, relaxation.relax_pressure
    )
end

mutable struct MeanResiduals
    resx_principle::Array{Float64, 1}
    resy_principle::Array{Float64, 1}
    resz_principle::Array{Float64, 1}
    resc_principle::Array{Float64, 1}
    resx_global::Array{Float64, 1}
    resy_global::Array{Float64, 1}
    resz_global::Array{Float64, 1}
    resc_global::Array{Float64, 1}
end

function MeanResiduals(nvcycles::Int64)::MeanResiduals
    return MeanResiduals(
        zeros(Float64, nvcycles),
        zeros(Float64, nvcycles),
        zeros(Float64, nvcycles),
        zeros(Float64, nvcycles),
        zeros(Float64, nvcycles),
        zeros(Float64, nvcycles),
        zeros(Float64, nvcycles),
        zeros(Float64, nvcycles)
    )
end

end # module