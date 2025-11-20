module MultigridVCycle

include("core/SmoothAndRestrict.jl")
include("core/SmoothAndProlongate.jl")

using Plots
using Printf
import EarthBox.ModelDataContainer.MultiGrids2dContainer: MultigridData
import ..MultigridDataManager: MultigridData2d, MultigridData3d
import ..MultigridDataManager.MultigridStructures: MeanResiduals, VcycleParameters, Counters
import ..LevelManager: LevelData, LevelData2d
import ..LevelManager: reset_solution_grids_to_zero!
import ..ViscosityScaling
import ..ResidualPlotting
import ..ResidualPlottingCM
import ..ArrayStats
import ..CalculateMeanResiduals: calculate_scaled_and_mean_residuals!
import ..Smoother: stokes_continuity3d_viscous_smoother!
import ..Smoother: stokes_continuity2d_viscous_smoother!
import .SmoothAndRestrict: smooth_and_restrict!
import .SmoothAndProlongate: smooth_and_prolongate!

const DEBUG = false
const DEBUG_WRITE_ARRAYS = false

function execute_multigrid_vcycles!(
    multigrid_data::MultigridData3d
)::Nothing
    use_multimulti = multigrid_data.vcycle.use_multimulti
    make_plots = multigrid_data.vcycle.make_plots
    nvcycles = multigrid_data.vcycle.nvcycles
    for ivcycle = 1:nvcycles
        update_ivcycle!(multigrid_data.counters, ivcycle)
        smooth_and_restrict!(multigrid_data)
        ΔRxL1, ΔRyL1, ΔRzL1, ΔRcL1 = smooth_and_prolongate!(multigrid_data)
        (
            ΔRxL1_scaled, ΔRyL1_scaled, ΔRzL1_scaled, ΔRcL1_scaled
        ) = calculate_mean_residuals!(multigrid_data, ΔRxL1, ΔRyL1, ΔRzL1, ΔRcL1)

        if make_plots
            plot_residuals(multigrid_data, ΔRxL1_scaled, ΔRyL1_scaled, ΔRzL1_scaled, ΔRcL1_scaled)
        end
        if !use_multimulti
            ViscosityScaling.rescale_viscosity_for_levels!(multigrid_data)
        else
            ViscosityScaling.rescale_viscosity_for_levels_using_multimulti!(multigrid_data)
        end
        if is_converged(multigrid_data)
            println(">>>> Multigrid iteration $ivcycle converged")
            break
        end
    end
    return nothing
end

function update_ivcycle!(counters::Counters, ivcycle::Int64)::Nothing
    counters.ivcycle = ivcycle
    return nothing
end

function execute_multigrid_vcycles!(
    multigrid_data::MultigridData2d
)::Nothing
    use_multimulti = multigrid_data.vcycle.use_multimulti
    make_plots = multigrid_data.vcycle.make_plots
    nvcycles = multigrid_data.vcycle.nvcycles
    for ivcycle = 1:nvcycles
        println(">>>> Multigrid iteration $ivcycle")
        update_ivcycle!(multigrid_data.counters, ivcycle)
        smooth_and_restrict!(multigrid_data)
        ΔRxL1, ΔRyL1, ΔRcL1 = smooth_and_prolongate!(multigrid_data)
        (
            ΔRxL1_scaled, ΔRyL1_scaled, ΔRcL1_scaled
        ) = calculate_mean_residuals!(multigrid_data, ΔRxL1, ΔRyL1, ΔRcL1)
        if make_plots
            plot_residuals(multigrid_data, ΔRxL1_scaled, ΔRyL1_scaled, ΔRcL1_scaled)
        end
        if !use_multimulti
            ViscosityScaling.rescale_viscosity_for_levels!(multigrid_data)
        else
            ViscosityScaling.rescale_viscosity_for_levels_using_multimulti!(multigrid_data)
        end
        if is_converged_2d(multigrid_data)
            println(">>>> Multigrid iteration $ivcycle converged")
            break
        end
    end
    return nothing
end

function calculate_mean_residuals!(
    multigrid_data::MultigridData3d,
    ΔRxL1::Array{Float64, 3},
    ΔRyL1::Array{Float64, 3},
    ΔRzL1::Array{Float64, 3},
    ΔRcL1::Array{Float64, 3}
)::Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}}
    ivcycle = multigrid_data.counters.ivcycle
    level_vector = multigrid_data.level_vector

    resx_principle = multigrid_data.mean_residuals.resx_principle
    resy_principle = multigrid_data.mean_residuals.resy_principle
    resz_principle = multigrid_data.mean_residuals.resz_principle
    resc_principle = multigrid_data.mean_residuals.resc_principle
    (
        ΔRxL1_scaled, ΔRyL1_scaled, ΔRzL1_scaled, ΔRcL1_scaled
    ) = calculate_scaled_and_mean_residuals!(
        ivcycle, multigrid_data.residual_scaling_factors.stokesscale, 
        multigrid_data.residual_scaling_factors.continscale, 
        level_vector[1], 
        ΔRxL1, ΔRyL1, ΔRzL1, ΔRcL1,
        resx_principle, resy_principle, resz_principle, resc_principle
        )

    return ΔRxL1_scaled, ΔRyL1_scaled, ΔRzL1_scaled, ΔRcL1_scaled
end

function calculate_mean_residuals!(
    multigrid_data::MultigridData2d,
    ΔRxL1::Array{Float64, 2},
    ΔRyL1::Array{Float64, 2},
    ΔRcL1::Array{Float64, 2}
)::Tuple{Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}}
    ivcycle = multigrid_data.counters.ivcycle
    level_vector = multigrid_data.level_vector

    resx_principle = multigrid_data.mean_residuals.resx_principle
    resy_principle = multigrid_data.mean_residuals.resy_principle
    resc_principle = multigrid_data.mean_residuals.resc_principle
    (
        ΔRxL1_scaled, ΔRyL1_scaled, ΔRcL1_scaled
    ) = calculate_scaled_and_mean_residuals!(
        ivcycle, multigrid_data.residual_scaling_factors.stokesscale, 
        multigrid_data.residual_scaling_factors.continscale, 
        level_vector[1], 
        ΔRxL1, ΔRyL1, ΔRcL1, 
        resx_principle, resy_principle, resc_principle 
        )

    return ΔRxL1_scaled, ΔRyL1_scaled, ΔRcL1_scaled
end

function plot_residuals(
    multigrid_data::MultigridData3d,
    ΔRxL1_scaled::Array{Float64, 3},
    ΔRyL1_scaled::Array{Float64, 3},
    ΔRzL1_scaled::Array{Float64, 3},
    ΔRcL1_scaled::Array{Float64, 3}
)::Nothing
    ivcycle = multigrid_data.counters.ivcycle
    level_vector = multigrid_data.level_vector
    nvcycles = multigrid_data.vcycle.nvcycles

    resx_principle = multigrid_data.mean_residuals.resx_principle
    resy_principle = multigrid_data.mean_residuals.resy_principle
    resz_principle = multigrid_data.mean_residuals.resz_principle
    resc_principle = multigrid_data.mean_residuals.resc_principle

    ResidualPlottingCM.plot_residual_surfaces(
        ΔRxL1_scaled, ΔRyL1_scaled, ΔRzL1_scaled, ΔRcL1_scaled, level_vector[1], ivcycle)
    ResidualPlottingCM.plot_residuals(resx_principle, resy_principle, resz_principle, resc_principle, nvcycles, ivcycle)
    return nothing
end

function plot_residuals(
    multigrid_data::MultigridData2d,
    ΔRxL1_scaled::Array{Float64, 2},
    ΔRyL1_scaled::Array{Float64, 2},
    ΔRcL1_scaled::Array{Float64, 2}
)::Nothing
    ivcycle = multigrid_data.counters.ivcycle
    level_vector = multigrid_data.level_vector
    nvcycles = multigrid_data.vcycle.nvcycles

    resx_principle = multigrid_data.mean_residuals.resx_principle
    resy_principle = multigrid_data.mean_residuals.resy_principle
    resc_principle = multigrid_data.mean_residuals.resc_principle

    ResidualPlottingCM.plot_residual_surfaces(
        ΔRxL1_scaled, ΔRyL1_scaled, ΔRcL1_scaled, level_vector[1], ivcycle)
    ResidualPlottingCM.plot_residuals_2d(
        resx_principle, resy_principle, resc_principle, nvcycles, ivcycle, ymin=-18.0, ymax=-1.0)
    return nothing
end

function is_converged(multigrid_data::MultigridData3d)::Bool
    ivcycle = multigrid_data.counters.ivcycle
    criterion = multigrid_data.vcycle.convergence_criterion
    if !multigrid_data.vcycle.use_multimulti
        max_principle_residual = calculate_max_principle_residual_3d(multigrid_data)
        nvcycles_viscosity_jump = multigrid_data.viscosity_scaling.nvcycles_viscosity_jump
        nviscosity_jumps = multigrid_data.viscosity_scaling.nviscosity_jumps
        nvcycles_to_max_viscosity = nvcycles_viscosity_jump * nviscosity_jumps
        bool_test = max_principle_residual < criterion && ivcycle > nvcycles_to_max_viscosity
        #if DEBUG
        println(">>>> ivcycle: $ivcycle, max_principle_residual: $max_principle_residual, criterion: $criterion, nvcycles_to_max_viscosity: $nvcycles_to_max_viscosity")
        #end
    else
        max_global_residual = calculate_max_global_residual_3d(multigrid_data)
        bool_test = max_global_residual < criterion
        #if DEBUG
        println(">>>> ivcycle: $ivcycle, max_global_residual: $max_global_residual, criterion: $criterion")
        #end
    end
    return bool_test
end

function is_converged_2d(multigrid_data::MultigridData2d)::Bool
    ivcycle = multigrid_data.counters.ivcycle
    criterion = multigrid_data.vcycle.convergence_criterion
    if !multigrid_data.vcycle.use_multimulti
        max_principle_residual = calculate_max_principle_residual_2d(multigrid_data)
        nvcycles_viscosity_jump = multigrid_data.viscosity_scaling.nvcycles_viscosity_jump
        nviscosity_jumps = multigrid_data.viscosity_scaling.nviscosity_jumps
        nvcycles_to_max_viscosity = nvcycles_viscosity_jump * nviscosity_jumps
        bool_test = max_principle_residual < criterion && ivcycle > nvcycles_to_max_viscosity
        #if DEBUG
        println(">>>> ivcycle: $ivcycle, max_principle_residual: $max_principle_residual, criterion: $criterion, nvcycles_to_max_viscosity: $nvcycles_to_max_viscosity, bool_test: $bool_test")
        #end
    else
        max_global_residual = calculate_max_global_residual_2d(multigrid_data)
        bool_test = max_global_residual < criterion
        #if DEBUG
        println(">>>> ivcycle: $ivcycle, max_global_residual: $max_global_residual, criterion: $criterion, bool_test: $bool_test")
        #end
    end
    return bool_test
end

function calculate_max_principle_residual_3d(multigrid_data::MultigridData3d)::Float64
    ivcycle = multigrid_data.counters.ivcycle
    resx_principle = multigrid_data.mean_residuals.resx_principle
    resy_principle = multigrid_data.mean_residuals.resy_principle
    resz_principle = multigrid_data.mean_residuals.resz_principle
    resc_principle = multigrid_data.mean_residuals.resc_principle
    if ivcycle > 2
        resx_principle = abs(10^resx_principle[ivcycle])
        resy_principle = abs(10^resy_principle[ivcycle])
        resz_principle = abs(10^resz_principle[ivcycle])
        resc_principle = abs(10^resc_principle[ivcycle])
        res_max = max(resx_principle, resy_principle, resz_principle, resc_principle)
        return res_max
    else
        return 1e38
    end
end

function calculate_max_principle_residual_2d(multigrid_data::MultigridData2d)::Float64
    ivcycle = multigrid_data.counters.ivcycle
    resx_principle = multigrid_data.mean_residuals.resx_principle
    resy_principle = multigrid_data.mean_residuals.resy_principle
    resc_principle = multigrid_data.mean_residuals.resc_principle
    if ivcycle > 2
        resx_principle = abs(10^resx_principle[ivcycle])
        resy_principle = abs(10^resy_principle[ivcycle])
        resc_principle = abs(10^resc_principle[ivcycle])
        res_max = max(resx_principle, resy_principle, resc_principle)
        return res_max
    else
        return 1e38
    end
end

function calculate_max_global_residual_3d(multigrid_data::MultigridData3d)::Float64
    iglobal_update = multigrid_data.counters.iglobal_update - 1
    resx_global = multigrid_data.mean_residuals.resx_global
    resy_global = multigrid_data.mean_residuals.resy_global
    resz_global = multigrid_data.mean_residuals.resz_global
    resc_global = multigrid_data.mean_residuals.resc_global
    if iglobal_update > 1
        resx_global = abs(10^resx_global[iglobal_update])
        resy_global = abs(10^resy_global[iglobal_update])
        resz_global = abs(10^resz_global[iglobal_update])
        resc_global = abs(10^resc_global[iglobal_update])
        res_max = max(resx_global, resy_global, resz_global, resc_global)
        return res_max
    else
        return 1e38
    end
end

function calculate_max_global_residual_2d(multigrid_data::MultigridData2d)::Float64
    iglobal_update = multigrid_data.counters.iglobal_update - 1
    resx_global = multigrid_data.mean_residuals.resx_global
    resy_global = multigrid_data.mean_residuals.resy_global
    resc_global = multigrid_data.mean_residuals.resc_global
    if iglobal_update > 1
        resx_global = abs(10^resx_global[iglobal_update])
        resy_global = abs(10^resy_global[iglobal_update])
        resc_global = abs(10^resc_global[iglobal_update])
        res_max = max(resx_global, resy_global, resc_global)
        return res_max
    else
        return 1e38
    end
end

end # module