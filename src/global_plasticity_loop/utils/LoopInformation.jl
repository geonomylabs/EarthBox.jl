module LoopInformation

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.EarthBoxDtypes: GlobalIterDictType
import Printf: @sprintf

function process_global_loop_info_and_update_dict!(
    model::ModelData,
    iglobal::Int,
    nglobal::Int,
    global_yield_error::Float64
)::Nothing
    main_time_loop = model.timestep.parameters.main_time_loop
    picard = model.stokes_continuity.parameters.picard
    conversion = model.conversion.parameters
    solution_norms = model.stokes_continuity.parameters.solution_norms
    residual_norms = model.stokes_continuity.parameters.residual_norms

    sec_per_myr = conversion.sec_per_Myr.value
    timestep = main_time_loop.timestep.value
    timesum = main_time_loop.timesum.value

    sec_per_myr = model.conversion.parameters.sec_per_Myr.value
    timesum_myr = timesum/sec_per_myr

    dvxy_rel_L2 = solution_norms.dvxy_rel_L2.value
    dvxy_rel_inf = solution_norms.dvxy_rel_inf.value
    dpr1_rel_L2 = solution_norms.dpr1_rel_L2.value
    dsoluv1_rel_L2 = solution_norms.dsoluv1_rel_L2.value

    resnl_rel_L2 = residual_norms.resnl_rel_L2.value

    tolerance = picard.tolerance_picard.value
    timestep_yr = timestep/sec_per_myr*1e6

    if dvxy_rel_L2 > tolerance
        status = "Not_Converged"
    else
        status = "Convergence_Successful"
    end
    if iglobal == nglobal && dvxy_rel_L2 > tolerance
        status = "Convergence_Failed_(iteration_limit_reached)"
    end

    print_info(
        "Stokes Iteration: " *
        "iglobal $iglobal : " *
        "dvxy_rel_L2 $(@sprintf("%.6f", dvxy_rel_L2)) : " *
        "global_yield_error (Pa) $(@sprintf("%.5E", global_yield_error)) : " *
        "Tolerance $(@sprintf("%.5f", tolerance)) : " *
        "Timestep_yr $timestep_yr : " *
        "Status -> $status",
        level=2
    )

    if !haskey(model.global_iter_dict, "tolerance")
        model.global_iter_dict["tolerance"] = tolerance
        model.global_iter_dict["time_step_yr"] = timestep_yr
        model.global_iter_dict["headers"] = [
            "iglobal",
            "dvxy_rel_L2",
            "dvxy_rel_inf",
            "resnl_rel_L2",
            "dpr1_rel_L2",
            "dsoluv1_rel_L2",
            "global_yield_error",
            "time_Myr",
            "status"
        ]
    end

    nkeys = length(keys(model.global_iter_dict))
    model.global_iter_dict[nkeys-3] = [
        iglobal,
        round(dvxy_rel_L2, digits=7),
        round(dvxy_rel_inf, digits=7),
        round(resnl_rel_L2, digits=7),
        round(dpr1_rel_L2, digits=7),
        round(dsoluv1_rel_L2, digits=7),
        global_yield_error,
        round(timesum_myr, digits=5),
        status
    ]

    return nothing
end

end # module 