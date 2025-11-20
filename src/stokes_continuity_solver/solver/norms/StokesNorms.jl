module StokesNorms

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox: MathTools
import LinearAlgebra: norm

function calculate_stokes_norms!(model::ModelData)
    @timeit_memit "Finished calculating Stokes norms" begin
      calculate_vxy_array!(model)
      calculate_solution_norms!(model)
      calculate_residual_norms!(model)
    end
end

""" Calculate vxy array.

# Updated Arrays from group model.stokes_continuity.arrays.velocity_solution
 - vxy.array::Vector{Float64}: ((ynum + 1)*xnum + (xnum + 1)*ynum))
    - Solution array with vx and vy solutions.
"""
function calculate_vxy_array!(model::ModelData)
    MathTools.compute_vxy_vec(
        model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array,
        model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array,
        model.grids.parameters.geometry.xnum.value,
        model.grids.parameters.geometry.ynum.value,
        model.stokes_continuity.arrays.velocity_solution.vxy.array
    )
end

""" Calculate Stokes solution norms.

# Updated parameters from group model.stokes_continuity.parameters.solution_norms
 - `dsoluv1_abs_inf.value::Float64`
   -  Infinity norm of the change in the Stokes solution.
 - `dsoluv1_rel_inf.value::Float64`
   -  Relative infinity norm of the change in the Stokes solution.
 - `dsoluv1_abs_L2.value::Float64`
   -  L2 of the change in the Stokes solution.
 - `dsoluv1_rel_L2.value::Float64`
   -  Relative L2 of the change in the Stokes solution.
 - `dvx1_abs_L2.value::Float64`
   -  L2 norm of the change in x-component velocity solution.
 - `dvx1_rel_L2.value::Float64`
   -  Relative L2 norm of the change in x-component velocity solution.
 - `dvy1_abs_L2.value::Float64`
   -  L2 norm of the change in y-component velocity solution.
 - `dvy1_rel_L2.value::Float64`
   -  Relative L2 norm of the change in y-component velocity solution.
 - `dpr1_abs_L2.value::Float64`
   -  L2 norm of the change in pressure solution.
 - `dpr1_rel_L2.value::Float64`
   -  Relative L2 norm of the change in pressure solution.
 - `dvxy_abs_inf.value::Float64`
   -  Infinity norm of the Stokes velocity solution array.
 - `dvxy_rel_inf.value::Float64`
   -  Relative infinity norm of the Stokes velocity solution array.
 - `dvxy_abs_L2.value::Float64`
   -  L2 norm of the Stokes velocity solution array.
 - `dvxy_rel_L2.value::Float64`
   -  Relative L2 norm of the Stokes velocity solution array.
"""
function calculate_solution_norms!(model::ModelData)
    (
        model.stokes_continuity.parameters.solution_norms.dsoluv1_abs_inf.value,
        model.stokes_continuity.parameters.solution_norms.dsoluv1_rel_inf.value
    ) = MathTools.get_inf_norms(
        model.stokes_continuity.arrays.stokes_solution.soluv1.array,
        model.stokes_continuity.arrays.stokes_solution.soluv1_old.array
    )

    (
        model.stokes_continuity.parameters.solution_norms.dsoluv1_abs_L2.value,
        model.stokes_continuity.parameters.solution_norms.dsoluv1_rel_L2.value
    ) = MathTools.get_L2_norms(
        model.stokes_continuity.arrays.stokes_solution.soluv1.array,
        model.stokes_continuity.arrays.stokes_solution.soluv1_old.array
    )

    (
        model.stokes_continuity.parameters.solution_norms.dvx1_abs_L2.value,
        model.stokes_continuity.parameters.solution_norms.dvx1_rel_L2.value
    ) = MathTools.get_L2_norms(
        model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array,
        model.stokes_continuity.arrays.staggered_grid_velocity.vx1_old.array
    )

    (
        model.stokes_continuity.parameters.solution_norms.dvy1_abs_L2.value,
        model.stokes_continuity.parameters.solution_norms.dvy1_rel_L2.value
    ) = MathTools.get_L2_norms(
        model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array,
        model.stokes_continuity.arrays.staggered_grid_velocity.vy1_old.array
    )

    (
        model.stokes_continuity.parameters.solution_norms.dpr1_abs_L2.value,
        model.stokes_continuity.parameters.solution_norms.dpr1_rel_L2.value
    ) = MathTools.get_L2_norms(
        model.stokes_continuity.arrays.pressure.pr1.array,
        model.stokes_continuity.arrays.pressure.pr1_old.array
    )

    (
        model.stokes_continuity.parameters.solution_norms.dvxy_abs_inf.value,
        model.stokes_continuity.parameters.solution_norms.dvxy_rel_inf.value
    ) = MathTools.get_inf_norms(
        model.stokes_continuity.arrays.velocity_solution.vxy.array,
        model.stokes_continuity.arrays.velocity_solution.vxy_old.array
    )

    (
        model.stokes_continuity.parameters.solution_norms.dvxy_abs_L2.value,
        model.stokes_continuity.parameters.solution_norms.dvxy_rel_L2.value
    ) = MathTools.get_L2_norms(
        model.stokes_continuity.arrays.velocity_solution.vxy.array,
        model.stokes_continuity.arrays.velocity_solution.vxy_old.array
    )
end

"""
    calculate_residual_norms(model)

Calculate residual Stokes norms.

# Updated parameters from group model.stokes_continuity.parameters.residual_norms:
 - `resnl_L2.value::Float64`
   -  L2 norm of the Stokes system residual array.
 - `resnl_L2_ini.value::Float64`
   -  Initial L2 norm of the Stokes system residual array.
 - `resnl_rel_L2.value::Float64`
   -  Relative L2 norm of the Stokes system residual array.
 - `resx_L2.value::Float64`
   -  L2 norm of the Stokes vx residual.
 - `resy_L2.value::Float64`
   -  L2 norm of the Stokes vy residual.
 - `resc_L2.value::Float64`
   -  L2 norm of the Stokes pressure residual.
"""
function calculate_residual_norms!(model::ModelData)
    resnl_L2 = norm(model.stokes_continuity.arrays.residuals.resnl.array)
    resnl_L2_ini = model.stokes_continuity.parameters.residual_norms.resnl_L2_ini.value

    if model.stokes_continuity.parameters.picard.iglobal.value == 0
        resnl_L2_ini = resnl_L2
    end

    resnl_rel_L2 = resnl_L2 != 0.0 ? resnl_L2/resnl_L2_ini : 1e38

    model.stokes_continuity.parameters.residual_norms.resnl_L2.value = resnl_L2
    model.stokes_continuity.parameters.residual_norms.resnl_L2_ini.value = resnl_L2_ini
    model.stokes_continuity.parameters.residual_norms.resnl_rel_L2.value = resnl_rel_L2
    model.stokes_continuity.parameters.residual_norms.resx_L2.value = 
        norm(model.stokes_continuity.arrays.residuals.resx1.array)
    model.stokes_continuity.parameters.residual_norms.resy_L2.value = 
        norm(model.stokes_continuity.arrays.residuals.resy1.array)
    model.stokes_continuity.parameters.residual_norms.resc_L2.value = 
        norm(model.stokes_continuity.arrays.residuals.resc1.array)
end

"""
    calculate_nonlinear_residual_norms(model)

Calculate nonlinear residual Stokes norms.

# Updated parameters from group model.stokes_continuity.parameters.residual_norms:
 - `resnlx_L2.value::Float64`
   -  Non-linear L2 norm of the Stokes vx residual.
 - `resnly_L2.value::Float64`
   -  Non-linear L2 norm of the Stokes vy residual.
 - `resnlc_L2.value::Float64`
   -  Non-linear L2 norm of the Stokes pressure residual.
 - `resnlx_L2_ini.value::Float64`
   -  Initial non-linear L2 norm of the Stokes vx residual.
 - `resnly_L2_ini.value::Float64`
   -  Initial non-linear L2 norm of the Stokes vy residual.
 - `resnlc_L2_ini.value::Float64`
   -  Initial non-linear L2 norm of the Stokes pressure residual.
 - `resnlx_rel_L2.value::Float64`
   -  Relative non-linear L2 norm of the Stokes vx residual.
 - `resnly_rel_L2.value::Float64`
   -  Relative non-linear L2 norm of the Stokes vy residual.
 - `resnlc_rel_L2.value::Float64`
   -  Relative non-linear L2 norm of the Stokes pressure residual.
"""
function calculate_nonlinear_residual_norms!(model::ModelData)
    residual_arrays = model.stokes_continuity.arrays.residuals
    residual_norms = model.stokes_continuity.parameters.residual_norms

    resnlx_L2 = norm(residual_arrays.resnlx.array)
    resnly_L2 = norm(residual_arrays.resnly.array)
    resnlc_L2 = norm(residual_arrays.resnlc.array)

    resnlx_L2_ini = residual_norms.resnlx_L2_ini.value
    resnly_L2_ini = residual_norms.resnly_L2_ini.value
    resnlc_L2_ini = residual_norms.resnlc_L2_ini.value

    if model.stokes_continuity.parameters.picard.iglobal.value == 0
        resnlx_L2_ini = resnlx_L2
        resnly_L2_ini = resnly_L2
        resnlc_L2_ini = resnlc_L2
    end

    resnlx_rel_L2 = resnlx_L2/resnlx_L2_ini
    resnly_rel_L2 = resnly_L2/resnly_L2_ini
    resnlc_rel_L2 = resnlc_L2/resnlc_L2_ini

    residual_norms.resnlx_L2.value = resnlx_L2
    residual_norms.resnly_L2.value = resnly_L2
    residual_norms.resnlc_L2.value = resnlc_L2
    residual_norms.resnlx_L2_ini.value = resnlx_L2_ini
    residual_norms.resnly_L2_ini.value = resnly_L2_ini
    residual_norms.resnlc_L2_ini.value = resnlc_L2_ini
    residual_norms.resnlx_rel_L2.value = resnlx_rel_L2
    residual_norms.resnly_rel_L2.value = resnly_rel_L2
    residual_norms.resnlc_rel_L2.value = resnlc_rel_L2
end

function print_all_norms(model::ModelData)
    output_level_2 = 1
    if output_level_2 == 1
        main_time_loop = model.timestep.parameters.main_time_loop
        picard = model.stokes_continuity.parameters.picard

        ntimestep = main_time_loop.ntimestep.value
        iglobal = picard.iglobal.value
        solution_norms = model.stokes_continuity.parameters.solution_norms
        residual_norms = model.stokes_continuity.parameters.residual_norms

        norm_obj_list = [
            solution_norms.dsoluv1_abs_inf,
            solution_norms.dsoluv1_rel_inf,
            solution_norms.dsoluv1_abs_L2,
            solution_norms.dsoluv1_rel_L2,
            solution_norms.dvx1_abs_L2,
            solution_norms.dvx1_rel_L2,
            solution_norms.dvy1_abs_L2,
            solution_norms.dvy1_rel_L2,
            solution_norms.dvxy_abs_inf,
            solution_norms.dvxy_rel_inf,
            solution_norms.dvxy_abs_L2,
            solution_norms.dvxy_rel_L2,
            solution_norms.dvxy_abs_L2,
            solution_norms.dvxy_rel_L2,
            residual_norms.resnl_L2,
            residual_norms.resnl_rel_L2,
            residual_norms.resx_L2,
            residual_norms.resy_L2,
            residual_norms.resc_L2,
            residual_norms.resnlx_L2,
            residual_norms.resnly_L2,
            residual_norms.resnlc_L2,
            residual_norms.resnlx_rel_L2,
            residual_norms.resnly_rel_L2,
            residual_norms.resnlc_rel_L2
        ]
        for norm_obj in norm_obj_list
            print_norm(ntimestep, iglobal, norm_obj)
        end
    end
end

function print_norm(ntimestep::Int, iglobal::Int, norm_obj)
    name = norm_obj.name
    norm = norm_obj.value
    println(
        "   NORMS : ntimestep : $ntimestep : iglobal : $iglobal : " *
        "$name : $norm"
    )
end

end # module 