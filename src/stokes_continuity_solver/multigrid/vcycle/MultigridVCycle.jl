module MultigridVCycle

using Printf
import EarthBox.ModelDataContainer.MultiGrids2dContainer: MultigridData
import ..MultigridDataManager: MultigridData2d, MultigridData3d
import ..MultigridDataManager.MultigridStructures: MeanResiduals, VcycleParameters, Counters
import ..LevelManager: LevelData, LevelData2d
import ..LevelManager: reset_solution_grids_to_zero!
import ..ViscosityScaling
import ..ResidualPlottingCM
import ..ArrayStats
import ..CalculateMeanResiduals: calculate_scaled_and_mean_residuals!
import ..CalculateMeanResiduals: accumulate_principle_residuals_3d!
import ..Smoother: stokes_continuity3d_viscous_smoother!
import ..Smoother: stokes_continuity2d_viscous_smoother!

const DEBUG = false
const DEBUG_WRITE_ARRAYS = false

"""When `EARTHBOX_MG_TIMING=1`, accumulate per-phase times over all V-cycles and print one summary."""
function mg_timing_enabled()::Bool
    return get(ENV, "EARTHBOX_MG_TIMING", "0") == "1"
end

"""When `EARTHBOX_MG_TIMING=1` and `EARTHBOX_MG_TIMING_DETAIL=1`, split restrict/prolong phases into smoother vs transfer."""
function mg_timing_detail_enabled()::Bool
    return mg_timing_enabled() && get(ENV, "EARTHBOX_MG_TIMING_DETAIL", "0") == "1"
end

const _mg_dt_sr_smooth = Ref(0.0)
const _mg_dt_sr_restrict = Ref(0.0)
const _mg_dt_sp_smooth = Ref(0.0)
const _mg_dt_sp_prolong = Ref(0.0)

function mg_timing_detail_reset!()::Nothing
    if mg_timing_detail_enabled()
        _mg_dt_sr_smooth[] = 0.0
        _mg_dt_sr_restrict[] = 0.0
        _mg_dt_sp_smooth[] = 0.0
        _mg_dt_sp_prolong[] = 0.0
    end
    return nothing
end

function mg_timing_detail_add_sr_smooth!(dt::Float64)::Nothing
    if mg_timing_detail_enabled()
        _mg_dt_sr_smooth[] += dt
    end
    return nothing
end

function mg_timing_detail_add_sr_restrict!(dt::Float64)::Nothing
    if mg_timing_detail_enabled()
        _mg_dt_sr_restrict[] += dt
    end
    return nothing
end

function mg_timing_detail_add_sp_smooth!(dt::Float64)::Nothing
    if mg_timing_detail_enabled()
        _mg_dt_sp_smooth[] += dt
    end
    return nothing
end

function mg_timing_detail_add_sp_prolong!(dt::Float64)::Nothing
    if mg_timing_detail_enabled()
        _mg_dt_sp_prolong[] += dt
    end
    return nothing
end

include("core/SmoothAndRestrict.jl")
include("core/SmoothAndProlongate.jl")

import .SmoothAndRestrict: smooth_and_restrict!
import .SmoothAndProlongate: smooth_and_prolongate!

function execute_multigrid_vcycles!(
    multigrid_data::MultigridData3d
)::Nothing
    use_multimulti = multigrid_data.vcycle.use_multimulti
    make_plots = multigrid_data.vcycle.make_plots
    nvcycles = multigrid_data.vcycle.nvcycles

    time_sr = 0.0
    time_sp = 0.0
    time_mr = 0.0
    time_vs = 0.0
    time_cv = 0.0
    mg_timing_detail_reset!()

    for ivcycle = 1:nvcycles
        update_ivcycle!(multigrid_data.counters, ivcycle)
        t0 = initialize_time()
        
        # Walk from fine to coarse levels, smoothing and restricting where restriction
        # involves interpolating residuals to coarser levels that are then used to update
        # the solution on the coarser levels via the smoother.
        smooth_and_restrict!(multigrid_data)
        time_sr, t0 = update_time(time_sr, t0)

        ΔRxL1, ΔRyL1, ΔRzL1, ΔRcL1 = smooth_and_prolongate!(multigrid_data)
        time_sp, t0 = update_time(time_sp, t0)

        if make_plots
            (
                ΔRxL1_scaled, ΔRyL1_scaled, ΔRzL1_scaled, ΔRcL1_scaled
            ) = calculate_mean_residuals!(multigrid_data, ΔRxL1, ΔRyL1, ΔRzL1, ΔRcL1)
            plot_residuals(multigrid_data, ΔRxL1_scaled, ΔRyL1_scaled, ΔRzL1_scaled, ΔRcL1_scaled)
        else
            calculate_mean_residuals_no_scaled_output!(multigrid_data, ΔRxL1, ΔRyL1, ΔRzL1, ΔRcL1)
        end
        time_mr, t0 = update_time(time_mr, t0)

        if !use_multimulti
            ViscosityScaling.rescale_viscosity_for_levels!(multigrid_data)
        else
            ViscosityScaling.rescale_viscosity_for_levels_using_multimulti!(multigrid_data)
        end
        time_vs, t0 = update_time(time_vs, t0)

        converged = is_converged(multigrid_data)
        time_cv, t0 = update_time(time_cv, t0)

        if converged
            println(">>>> Multigrid iteration $ivcycle converged")
            break
        end
    end

    print_timing_summary(time_sr, time_sp, time_mr, time_vs, time_cv)
    
    return nothing
end

function initialize_time()::Float64
    if mg_timing_enabled()
        return time()
    else
        return 0.0
    end
end

function update_time(time_sum::Float64, t0::Float64)::Tuple{Float64, Float64}
    if mg_timing_enabled()
        time_sum += time() - t0
        t0 = time()
        return time_sum, t0
    else
        return time_sum, t0
    end
end

function print_timing_summary(  
    time_sr::Float64, 
    time_sp::Float64, 
    time_mr::Float64,
    time_vs::Float64,
    time_cv::Float64,
)::Nothing
    if mg_timing_enabled()
        tot = time_sr + time_sp + time_mr + time_vs + time_cv
        println(">> MG timing (EARTHBOX_MG_TIMING=1): total tracked=$(tot)s")
        println(">>   smooth_and_restrict: $(time_sr)s ($(100 * time_sr / tot)% )")
        println(">>   smooth_and_prolongate: $(time_sp)s ($(100 * time_sp / tot)% )")
        println(">>   mean_residuals: $(time_mr)s ($(100 * time_mr / tot)% )")
        println(">>   viscosity_rescale: $(time_vs)s ($(100 * time_vs / tot)% )")
        println(">>   convergence_check: $(time_cv)s ($(100 * time_cv / tot)% )")
        if mg_timing_detail_enabled()
            inner = _mg_dt_sr_smooth[] + _mg_dt_sr_restrict[] +
                _mg_dt_sp_smooth[] + _mg_dt_sp_prolong[]
            pct(x) = inner > 0 ? 100 * x / inner : 0.0
            println(">> MG detail (EARTHBOX_MG_TIMING_DETAIL=1), summed over V-cycles:")
            println(
                ">>   SR smoother+residual: $(_mg_dt_sr_smooth[])s ($(pct(_mg_dt_sr_smooth[]))% of detail)",
            )
            println(
                ">>   SR restriction: $(_mg_dt_sr_restrict[])s ($(pct(_mg_dt_sr_restrict[]))% of detail)",
            )
            println(
                ">>   SP smoother+residual: $(_mg_dt_sp_smooth[])s ($(pct(_mg_dt_sp_smooth[]))% of detail)",
            )
            println(
                ">>   SP prolongation+correction: $(_mg_dt_sp_prolong[])s ($(pct(_mg_dt_sp_prolong[]))% of detail)",
            )
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
        tt1 = time()
        println("--------------------------------")
        println("> Multigrid iteration $ivcycle")
        println("--------------------------------")

        t1 = time()
        println(">> Updating ivcycle...")
        update_ivcycle!(multigrid_data.counters, ivcycle)
        t2 = time()
        println(">> Time taken to update ivcycle: $(t2-t1) seconds")

        t1 = time()
        println(">> Smoothing and restricting...")
        smooth_and_restrict!(multigrid_data)
        t2 = time()
        println(">> Time taken to smooth and restrict: $(t2-t1) seconds")

        t1 = time()
        println(">> Smoothing and prolongating...")
        ΔRxL1, ΔRyL1, ΔRcL1 = smooth_and_prolongate!(multigrid_data)
        t2 = time()
        println(">> Time taken to smooth and prolongate: $(t2-t1) seconds")

        t1 = time()
        println(">> Calculating mean residuals...")
        (
            ΔRxL1_scaled, ΔRyL1_scaled, ΔRcL1_scaled
        ) = calculate_mean_residuals!(multigrid_data, ΔRxL1, ΔRyL1, ΔRcL1)
        t2 = time()
        println(">> Time taken to calculate mean residuals: $(t2-t1) seconds")

        t1 = time()
        println(">> Plotting residuals...")
        if make_plots
            plot_residuals(multigrid_data, ΔRxL1_scaled, ΔRyL1_scaled, ΔRcL1_scaled)
        end
        t2 = time()
        println(">> Time taken to plot residuals: $(t2-t1) seconds")

        t1 = time()
        println(">> Rescaling viscosity...")
        if !use_multimulti
            ViscosityScaling.rescale_viscosity_for_levels!(multigrid_data)
        else
            ViscosityScaling.rescale_viscosity_for_levels_using_multimulti!(multigrid_data)
        end
        t2 = time()
        println(">> Time taken to rescale viscosity: $(t2-t1) seconds")

        t1 = time()
        println(">> Checking convergence...")
        if is_converged_2d(multigrid_data)
            println("> Multigrid iteration $ivcycle converged")
            break
        end
        t2 = time()
        println(">> Time taken to check convergence: $(t2-t1) seconds")

        tt2 = time()
        println(">>> Total time taken for iteration $ivcycle: $(tt2-tt1) seconds")

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

function calculate_mean_residuals_no_scaled_output!(
    multigrid_data::MultigridData3d,
    ΔRxL1::Array{Float64, 3},
    ΔRyL1::Array{Float64, 3},
    ΔRzL1::Array{Float64, 3},
    ΔRcL1::Array{Float64, 3}
)::Nothing
    ivcycle = multigrid_data.counters.ivcycle
    level_vector = multigrid_data.level_vector
    accumulate_principle_residuals_3d!(
        ivcycle,
        multigrid_data.residual_scaling_factors.stokesscale,
        multigrid_data.residual_scaling_factors.continscale,
        level_vector[1],
        ΔRxL1,
        ΔRyL1,
        ΔRzL1,
        ΔRcL1,
        multigrid_data.mean_residuals.resx_principle,
        multigrid_data.mean_residuals.resy_principle,
        multigrid_data.mean_residuals.resz_principle,
        multigrid_data.mean_residuals.resc_principle,
    )
    return nothing
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
        println(">>>> ivcycle: $ivcycle, max_principle_residual: $max_principle_residual, criterion: $criterion, nvcycles_to_max_viscosity: $nvcycles_to_max_viscosity")
    else
        # Multimulti: use global residual only (same as 2D multimulti). Principle RMS is updated every
        # V-cycle and can show sawtooth spikes that do not indicate sustained convergence; global
        # residuals are recorded at global viscosity-update steps and reflect continuation progress.
        max_global_residual = calculate_max_global_residual_3d(multigrid_data)
        bool_test = max_global_residual < criterion
        println(
            ">>>> ivcycle: $ivcycle, max_global_residual: $max_global_residual, " *
            "criterion: $criterion, bool_test: $bool_test",
        )
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