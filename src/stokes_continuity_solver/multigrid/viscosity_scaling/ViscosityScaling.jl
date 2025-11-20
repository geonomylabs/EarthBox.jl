module ViscosityScaling

import ..MultigridDataManager: MultigridData3d, MultigridData2d
import ..MultigridDataManager.MultigridStructures: Counters, ViscosityScalingParameters
import ..LevelManager: LevelData, LevelData2d
import ..ResidualPlotting
import ..ResidualPlottingCM
import ..ArrayStats
import ..Smoother: update_rhs_parts_on_level1_using_global_level0_residuals
import ..CalculateMeanResiduals: calculate_scaled_and_mean_residuals!

function rescale_viscosity!(
    multigrid_data::MultigridData3d
)::Nothing
    viscosity_scaling = multigrid_data.viscosity_scaling
    level_vector = multigrid_data.level_vector
    viscosity_max_o = viscosity_scaling.viscosity_max_o
    viscosity_min_o = viscosity_scaling.viscosity_min_o
    viscosity_max_cur = viscosity_scaling.viscosity_max_cur
    viscosity_min_cur = viscosity_scaling.viscosity_min_cur
    if viscosity_max_cur < viscosity_max_o || viscosity_min_cur > viscosity_min_o
        scale_viscosity_for_all_levels!(level_vector, viscosity_scaling)
    end
    return nothing
end

function rescale_viscosity!(
    multigrid_data::MultigridData2d
)::Nothing
    viscosity_scaling = multigrid_data.viscosity_scaling
    level_vector = multigrid_data.level_vector
    viscosity_max_o = viscosity_scaling.viscosity_max_o
    viscosity_min_o = viscosity_scaling.viscosity_min_o
    viscosity_max_cur = viscosity_scaling.viscosity_max_cur
    viscosity_min_cur = viscosity_scaling.viscosity_min_cur
    if viscosity_max_cur < viscosity_max_o || viscosity_min_cur > viscosity_min_o
        scale_viscosity_for_all_levels!(level_vector, viscosity_scaling)
    end
    return nothing
end

function scale_viscosity_for_all_levels!(
    level_vector::Vector{LevelData},
    viscosity_scaling::ViscosityScalingParameters
)::Nothing
    viscosity_max_o = viscosity_scaling.viscosity_max_o
    viscosity_min_o = viscosity_scaling.viscosity_min_o
    viscosity_min_cur = viscosity_scaling.viscosity_min_cur
    viscosity_max_cur = viscosity_scaling.viscosity_max_cur
    eta_coeff = log(viscosity_max_cur / viscosity_min_cur) / log(viscosity_max_o / viscosity_min_o)
    levelnum = length(level_vector)
    for n = 1:levelnum
        level_vector[n].etaxy.array .= viscosity_min_cur * exp.(log.(level_vector[n].etaxyo.array ./ viscosity_min_o) * eta_coeff)
        level_vector[n].etaxz.array .= viscosity_min_cur * exp.(log.(level_vector[n].etaxzo.array ./ viscosity_min_o) * eta_coeff)
        level_vector[n].etayz.array .= viscosity_min_cur * exp.(log.(level_vector[n].etayzo.array ./ viscosity_min_o) * eta_coeff)
        level_vector[n].etan.array .= viscosity_min_cur * exp.(log.(level_vector[n].etano.array ./ viscosity_min_o) * eta_coeff)
    end
    return nothing
end

function scale_viscosity_for_all_levels!(
    level_vector::Vector{LevelData2d},
    viscosity_scaling::ViscosityScalingParameters
)::Nothing
    viscosity_max_o = viscosity_scaling.viscosity_max_o
    viscosity_min_o = viscosity_scaling.viscosity_min_o
    viscosity_min_cur = viscosity_scaling.viscosity_min_cur
    viscosity_max_cur = viscosity_scaling.viscosity_max_cur
    eta_coeff = log(viscosity_max_cur / viscosity_min_cur) / log(viscosity_max_o / viscosity_min_o)
    levelnum = length(level_vector)
    for n = 1:levelnum
        level_vector[n].etas.array .= viscosity_min_cur * exp.(log.(level_vector[n].etaso.array ./ viscosity_min_o) * eta_coeff)
        level_vector[n].etan.array .= viscosity_min_cur * exp.(log.(level_vector[n].etano.array ./ viscosity_min_o) * eta_coeff)
    end
    return nothing
end

function rescale_viscosity_for_levels!(
    multigrid_data::MultigridData3d
)::Nothing
    ivcycle = multigrid_data.counters.ivcycle
    level_vector = multigrid_data.level_vector

    viscosity_max_cur = multigrid_data.viscosity_scaling.viscosity_max_cur
    viscosity_min_cur = multigrid_data.viscosity_scaling.viscosity_min_cur
    viscosity_max_o = multigrid_data.viscosity_scaling.viscosity_max_o
    viscosity_min_o = multigrid_data.viscosity_scaling.viscosity_min_o
    nvcycles_viscosity_jump = multigrid_data.viscosity_scaling.nvcycles_viscosity_jump
    resmax = multigrid_data.viscosity_scaling.resmax

    resx_principle = multigrid_data.mean_residuals.resx_principle
    # Renormalizing viscosity
    # Add cycle counter
    iviscosity_jump = increment_iviscosity_jump!(multigrid_data.counters)
    # In Gerya's example resx_principle is used for all residuals. This seems odd.
    # Consider using the following formula for restot:
    # restot = (10^resx_principle[ivcycle] + 10^resy_principle[ivcycle] + 10^resz_principle[ivcycle] + 10^resc_principle[ivcycle])
    # But the original formula for the code will be used for benchmarking
    restot = (10^resx_principle[ivcycle] + 10^resx_principle[ivcycle] + 10^resx_principle[ivcycle] + 10^resx_principle[ivcycle])
    if (iviscosity_jump > nvcycles_viscosity_jump || restot < resmax) && (viscosity_max_cur < viscosity_max_o || viscosity_min_cur > viscosity_min_o)
        reset_iviscosity_jump_to_one!(multigrid_data.counters)
        update_viscosity_min_cur!(multigrid_data.viscosity_scaling)
        update_viscosity_max_cur!(multigrid_data.viscosity_scaling)
        scale_viscosity_for_all_levels!(level_vector, multigrid_data.viscosity_scaling)
    end
    return nothing
end

function rescale_viscosity_for_levels!(
    multigrid_data::MultigridData2d,
)::Nothing
    level_vector = multigrid_data.level_vector

    viscosity_max_o = multigrid_data.viscosity_scaling.viscosity_max_o
    viscosity_min_o = multigrid_data.viscosity_scaling.viscosity_min_o
    nvcycles_viscosity_jump = multigrid_data.viscosity_scaling.nvcycles_viscosity_jump
    viscosity_max_cur = multigrid_data.viscosity_scaling.viscosity_max_cur
    viscosity_min_cur = multigrid_data.viscosity_scaling.viscosity_min_cur

    iviscosity_jump = increment_iviscosity_jump!(multigrid_data.counters)
    if (iviscosity_jump > nvcycles_viscosity_jump) && (viscosity_max_cur < viscosity_max_o || viscosity_min_cur > viscosity_min_o)
        reset_iviscosity_jump_to_one!(multigrid_data.counters)
        update_viscosity_min_cur!(multigrid_data.viscosity_scaling)
        update_viscosity_max_cur!(multigrid_data.viscosity_scaling)
        scale_viscosity_for_all_levels!(level_vector, multigrid_data.viscosity_scaling)
    end
    return nothing
end

function rescale_viscosity_for_levels_using_multimulti!(
    multigrid_data::MultigridData3d
)::Nothing
    level0 = multigrid_data.level0
    level_vector = multigrid_data.level_vector
    make_plots = multigrid_data.vcycle.make_plots

    nvcycles_viscosity_jump = multigrid_data.viscosity_scaling.nvcycles_viscosity_jump
    iviscosity_jump = increment_iviscosity_jump!(multigrid_data.counters)
    if iviscosity_jump > nvcycles_viscosity_jump
        reset_iviscosity_jump_to_one!(multigrid_data.counters)
        update_viscosity_max_cur_and_iviscosity_max_overstep!(
            multigrid_data.viscosity_scaling, multigrid_data.counters)
        iviscosity_max_overstep = multigrid_data.counters.iviscosity_max_overstep
        if iviscosity_max_overstep > 1
            reset_viscosity_max_to_starting_value!(multigrid_data.viscosity_scaling)
            reset_iviscosity_max_overstep_to_zero!(multigrid_data.counters)
            update_level0_solution_using_corrections_from_level1!(level0, level_vector[1])
            update_level0_viscosity_using_original!(level0, level_vector)
            (
                ΔRx_L0, ΔRy_L0, ΔRz_L0, ΔRc_L0
            ) = update_rhs_parts_on_level1_using_global_level0_residuals(multigrid_data)
            (
                ΔRx_L0_scaled, ΔRy_L0_scaled, ΔRz_L0_scaled, ΔRc_L0_scaled
            ) = calculate_global_mean_residuals!(multigrid_data, ΔRx_L0, ΔRy_L0, ΔRz_L0, ΔRc_L0)
            if make_plots
                plot_residuals(
                    multigrid_data, ΔRx_L0_scaled, ΔRy_L0_scaled, ΔRz_L0_scaled, ΔRc_L0_scaled)
            end
            set_level1_solution_to_zero!(level_vector[1])
            update_iglobal_update!(multigrid_data.counters)
        end
        scale_viscosity_for_all_levels!(level_vector, multigrid_data.viscosity_scaling)
    end
    return nothing
end

function rescale_viscosity_for_levels_using_multimulti!(
    multigrid_data::MultigridData2d
)::Nothing
    level0 = multigrid_data.level0
    level_vector = multigrid_data.level_vector
    make_plots = multigrid_data.vcycle.make_plots

    nvcycles_viscosity_jump = multigrid_data.viscosity_scaling.nvcycles_viscosity_jump
    iviscosity_jump = increment_iviscosity_jump!(multigrid_data.counters)
    if iviscosity_jump > nvcycles_viscosity_jump
        reset_iviscosity_jump_to_one!(multigrid_data.counters)
        update_viscosity_max_cur_and_iviscosity_max_overstep!(
            multigrid_data.viscosity_scaling, multigrid_data.counters)
        iviscosity_max_overstep = multigrid_data.counters.iviscosity_max_overstep
        if iviscosity_max_overstep > 1
            reset_viscosity_max_to_starting_value!(multigrid_data.viscosity_scaling)
            reset_iviscosity_max_overstep_to_zero!(multigrid_data.counters)
            update_level0_solution_using_corrections_from_level1!(level0, level_vector[1])
            update_level0_viscosity_using_original!(level0, level_vector)
            (
                ΔRx_L0, ΔRy_L0, ΔRc_L0
            ) = update_rhs_parts_on_level1_using_global_level0_residuals(multigrid_data)
            (
                ΔRx_L0_scaled, ΔRy_L0_scaled, ΔRc_L0_scaled
            ) = calculate_global_mean_residuals!(multigrid_data, ΔRx_L0, ΔRy_L0, ΔRc_L0)
            if make_plots
                plot_residuals(
                    multigrid_data, ΔRx_L0_scaled, ΔRy_L0_scaled, ΔRc_L0_scaled)
            end
            set_level1_solution_to_zero!(level_vector[1])
            update_iglobal_update!(multigrid_data.counters)
        end
        scale_viscosity_for_all_levels!(level_vector, multigrid_data.viscosity_scaling)
    end
    return nothing
end

function update_level0_solution_using_corrections_from_level1!(
    level0::LevelData2d,
    level1::LevelData2d
)::Nothing
    level0.vx.array .+= level1.vx.array
    level0.vy.array .+= level1.vy.array
    level0.pr.array .+= level1.pr.array
    return nothing
end

function update_level0_solution_using_corrections_from_level1!(
    level0::LevelData,
    level1::LevelData
)::Nothing
    level0.vx.array .+= level1.vx.array
    level0.vy.array .+= level1.vy.array
    level0.vz.array .+= level1.vz.array
    level0.pr.array .+= level1.pr.array
    return nothing
end

function update_level0_viscosity_using_original!(
    level0::LevelData2d,
    level_vector::Vector{LevelData2d}
)::Nothing
    level0.etan.array .= level_vector[1].etano.array
    level0.etas.array .= level_vector[1].etaso.array
    return nothing
end

function update_level0_viscosity_using_original!(
    level0::LevelData,
    level_vector::Vector{LevelData}
)::Nothing
    level0.etan.array .= level_vector[1].etano.array
    level0.etaxy.array .= level_vector[1].etaxyo.array
    level0.etaxz.array .= level_vector[1].etaxzo.array
    level0.etayz.array .= level_vector[1].etayzo.array
    return nothing
end

function set_level1_solution_to_zero!(
    level1::LevelData2d
)::Nothing
    level1.vx.array .= 0.0
    level1.vy.array .= 0.0
    level1.pr.array .= 0.0
    return nothing
end

function set_level1_solution_to_zero!(
    level1::LevelData
)::Nothing
    level1.vx.array .= 0.0
    level1.vy.array .= 0.0
    level1.vz.array .= 0.0
    level1.pr.array .= 0.0
    return nothing
end

function increment_iviscosity_jump!(
    counters::Counters
)::Int64
    counters.iviscosity_jump += 1
    return counters.iviscosity_jump
end

function reset_iviscosity_jump_to_one!(
    counters::Counters
)::Nothing
    counters.iviscosity_jump = 1
    return nothing
end

function update_viscosity_min_cur!(
    viscosity_scaling::ViscosityScalingParameters
)::Nothing
    viscosity_min_o = viscosity_scaling.viscosity_min_o
    viscosity_min_cur = viscosity_scaling.viscosity_min_cur
    viscosity_min_factor = viscosity_scaling.viscosity_min_factor

    viscosity_min_cur = viscosity_min_cur * viscosity_min_factor
    if viscosity_min_cur < viscosity_min_o
        viscosity_min_cur = viscosity_min_o
    end
    viscosity_scaling.viscosity_min_cur = viscosity_min_cur
    return nothing
end

function update_viscosity_max_cur!(
    viscosity_scaling::ViscosityScalingParameters
)::Nothing
    viscosity_max_o = viscosity_scaling.viscosity_max_o
    viscosity_max_cur = viscosity_scaling.viscosity_max_cur
    viscosity_max_factor = viscosity_scaling.viscosity_max_factor
    
    viscosity_max_cur = viscosity_max_cur * viscosity_max_factor
    if viscosity_max_cur > viscosity_max_o
        viscosity_max_cur = viscosity_max_o
    end
    viscosity_scaling.viscosity_max_cur = viscosity_max_cur
    return nothing
end

function update_viscosity_max_cur_and_iviscosity_max_overstep!(
    viscosity_scaling::ViscosityScalingParameters,
    counters::Counters,
)::Nothing
    viscosity_max_o = viscosity_scaling.viscosity_max_o
    viscosity_max_cur = viscosity_scaling.viscosity_max_cur
    viscosity_max_factor = viscosity_scaling.viscosity_max_factor
    
    viscosity_max_cur = viscosity_max_cur * viscosity_max_factor
    if viscosity_max_cur >= viscosity_max_o
        counters.iviscosity_max_overstep += 1
    end
    if viscosity_max_cur > viscosity_max_o
        viscosity_max_cur = viscosity_max_o
    end
    viscosity_scaling.viscosity_max_cur = viscosity_max_cur
    return nothing
end

function reset_viscosity_max_to_starting_value!(
    viscosity_scaling::ViscosityScalingParameters
)::Nothing
    viscosity_max_start = viscosity_scaling.viscosity_max_start
    viscosity_scaling.viscosity_max_cur = viscosity_max_start
    return nothing
end

function reset_iviscosity_max_overstep_to_zero!(
    counters::Counters
)::Nothing
    counters.iviscosity_max_overstep = 0
    return nothing
end

function update_iglobal_update!(
    counters::Counters,
)::Nothing
    counters.iglobal_update += 1
    return nothing
end

function calculate_global_mean_residuals!(
    multigrid_data::MultigridData3d,
    ΔRx_L0::Array{Float64, 3},
    ΔRy_L0::Array{Float64, 3},
    ΔRz_L0::Array{Float64, 3},
    ΔRc_L0::Array{Float64, 3}
)::Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}}
    iglobal_update = multigrid_data.counters.iglobal_update
    level_vector = multigrid_data.level_vector

    resx_global = multigrid_data.mean_residuals.resx_global
    resy_global = multigrid_data.mean_residuals.resy_global
    resz_global = multigrid_data.mean_residuals.resz_global
    resc_global = multigrid_data.mean_residuals.resc_global
    (
        ΔRx_L0_scaled, ΔRy_L0_scaled, ΔRz_L0_scaled, ΔRc_L0_scaled
    ) = calculate_scaled_and_mean_residuals!(
        iglobal_update, multigrid_data.residual_scaling_factors.stokesscale, 
        multigrid_data.residual_scaling_factors.continscale, 
        level_vector[1], 
        ΔRx_L0, ΔRy_L0, ΔRz_L0, ΔRc_L0,
        resx_global, resy_global, resz_global, resc_global
        )

    return ΔRx_L0_scaled, ΔRy_L0_scaled, ΔRz_L0_scaled, ΔRc_L0_scaled
end

function calculate_global_mean_residuals!(
    multigrid_data::MultigridData2d,
    ΔRx_L0::Array{Float64, 2},
    ΔRy_L0::Array{Float64, 2},
    ΔRc_L0::Array{Float64, 2}
)::Tuple{Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}}
    iglobal_update = multigrid_data.counters.iglobal_update
    level_vector = multigrid_data.level_vector

    resx_global = multigrid_data.mean_residuals.resx_global
    resy_global = multigrid_data.mean_residuals.resy_global
    resc_global = multigrid_data.mean_residuals.resc_global

    (
        ΔRx_L0_scaled, ΔRy_L0_scaled, ΔRc_L0_scaled
    ) = calculate_scaled_and_mean_residuals!(
        iglobal_update, multigrid_data.residual_scaling_factors.stokesscale, 
        multigrid_data.residual_scaling_factors.continscale, 
        level_vector[1], 
        ΔRx_L0, ΔRy_L0, ΔRc_L0,
        resx_global, resy_global, resc_global
        )

    return ΔRx_L0_scaled, ΔRy_L0_scaled, ΔRc_L0_scaled
end

function plot_residuals(
    multigrid_data::MultigridData3d,
    ΔRx_L0_scaled::Array{Float64, 3},
    ΔRy_L0_scaled::Array{Float64, 3},
    ΔRz_L0_scaled::Array{Float64, 3},
    ΔRc_L0_scaled::Array{Float64, 3}
)::Nothing
    level0 = multigrid_data.level0
    iglobal_update = multigrid_data.counters.iglobal_update
    nvcycles = multigrid_data.vcycle.nvcycles

    resx_global = multigrid_data.mean_residuals.resx_global
    resy_global = multigrid_data.mean_residuals.resy_global
    resz_global = multigrid_data.mean_residuals.resz_global
    resc_global = multigrid_data.mean_residuals.resc_global

    ResidualPlottingCM.plot_residual_surfaces(
        ΔRx_L0_scaled, ΔRy_L0_scaled, ΔRz_L0_scaled, ΔRc_L0_scaled, level0, 
        iglobal_update, msg="GlobalLevel0"
        )
    ResidualPlottingCM.plot_residuals(
        resx_global, resy_global, resz_global, resc_global, nvcycles, 
        iglobal_update, msg="GlobalLevel0", xlabel="Global Update Step"
        )
    return nothing
end

function plot_residuals(
    multigrid_data::MultigridData2d,
    ΔRx_L0_scaled::Array{Float64, 2},
    ΔRy_L0_scaled::Array{Float64, 2},
    ΔRc_L0_scaled::Array{Float64, 2}
)::Nothing
    level0 = multigrid_data.level0
    iglobal_update = multigrid_data.counters.iglobal_update
    nvcycles = multigrid_data.vcycle.nvcycles

    resx_global = multigrid_data.mean_residuals.resx_global
    resy_global = multigrid_data.mean_residuals.resy_global
    resc_global = multigrid_data.mean_residuals.resc_global

    ResidualPlottingCM.plot_residual_surfaces(
        ΔRx_L0_scaled, ΔRy_L0_scaled, ΔRc_L0_scaled, level0, 
        iglobal_update, msg="GlobalLevel0"
        )
    ResidualPlottingCM.plot_residuals_2d(
        resx_global, resy_global, resc_global, nvcycles, 
        iglobal_update, msg="GlobalLevel0", xlabel="Global Update Step"
        )
    return nothing
end

end