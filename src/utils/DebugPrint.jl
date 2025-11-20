module DebugPrint

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Arrays: ArrayUtils
import Printf: @printf
import Printf: @sprintf

const i_check = 25
const j_check = 1
const iglobal_check = 108795

function check_RHS(iglobal::Int64, loc::String)
    iglobal_check2 = 108795
    if iglobal == iglobal_check2
        println("**=> DEBUG check RHS: iglobal = $iglobal in $loc")
    end
    return nothing
end

function check_RHS_value(model::ModelData, loc::String)::Nothing
    RHS = model.stokes_continuity.arrays.rhs.RHS.array
    println("Check RHS at $iglobal_check : ", RHS[iglobal_check], " in $loc")
    return nothing
end

function check_RC_and_pscale(
    iglobal::Int64,
    r_value::Float64,
    rc_value::Float64, 
    pscale::Float64, 
    loc::String
)::Nothing
    if iglobal == iglobal_check
        println("**=> DEBUG check RC and pscale: iglobal = $iglobal in $loc")
        println("        > r_value = $r_value, rc_value = $rc_value, pscale = $pscale")
    end
    return nothing
end

function check_viscoelastic_normal_stress_and_viscosity_terms(
    i::Int64, 
    j::Int64,
    timestep::Float64, 
    etan0::Matrix{Float64},
    etan1::Matrix{Float64},
    mun1::Matrix{Float64},
)::Nothing
    if i == i_check && j == j_check
        println("=> DEBUG check viscoelastic terms calc: i = $i, j = $j, timestep = $timestep, etan0 = $(etan0[i, j]), etan1 = $(etan1[i, j]), mun1 = $(mun1[i, j])")
    end
    return nothing
end

function print_rhs_terms_along_middle(
    model::ModelData, 
    msg::String,
    print_row::Bool,
    print_column::Bool
)::Nothing
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array
    gridx_pr = model.grids.arrays.pressure.gridx_pr.array
    gridy_pr = model.grids.arrays.pressure.gridy_pr.array
    # Vx-grid: y coordinates of vx grid nodes; note that x-coordinates are gridx_b
    gridy_vx = model.grids.arrays.staggered_vx.gridy_vx.array 
    # Vy-grid: x coordinates of vy grid nodes; note that y coordinates are gridy_b
    gridx_vy = model.grids.arrays.staggered_vy.gridx_vy.array 

    RC1 = model.stokes_continuity.arrays.rhs.RC1.array
    RX1 = model.stokes_continuity.arrays.rhs.RX1.array
    RY1 = model.stokes_continuity.arrays.rhs.RY1.array

    println("DEBUG: printing rhs terms:", msg)

    ArrayUtils.print_min_max("RC1", RC1)
    ArrayUtils.print_min_max("RX1", RX1)
    ArrayUtils.print_min_max("RY1", RY1)

    if print_row
        i_middle_row = size(RC1, 1) ÷ 2 + 1
        println("==> RC1 along row $i_middle_row at y = $(gridy_pr[i_middle_row])")
        n = size(RC1, 2)
        for j in 1:n
            println("j = ", j, " x = ", gridx_pr[j], " RC1 = ", RC1[i_middle_row, j])
        end
        println()
        i_middle_row = size(RX1, 1) ÷ 2 + 1
        println("==> RX1 along row $i_middle_row at y = $(gridy_vx[i_middle_row])")
        n = size(RX1, 2)
        for j in 1:n
            println("j = ", j, " x = ", gridx_b[j], " RX1 = ", RX1[i_middle_row, j])
        end
        println()
        i_middle_row = size(RY1, 1) ÷ 2 + 1
        println("==> RY1 along row $i_middle_row at y = $(gridy_b[i_middle_row])")
        n = size(RY1, 2)
        for j in 1:n
            println("j = ", j, " x = ", gridx_vy[j], " RY1 = ", RY1[i_middle_row, j])
        end
        println()
    end

    if print_column
        j_middle_column = size(RC1, 2) ÷ 2 + 1
        println("==> RC1 along column $j_middle_column at x = $(gridx_pr[j_middle_column])")
        n = size(RC1, 1)
        for i in 1:n
            println("i = ", i, " y = ", gridy_pr[i], " RC1 = ", RC1[i, j_middle_column])
        end
        println()
        j_middle_column = size(RX1, 2) ÷ 2 + 1
        println("==> RX1 along column $j_middle_column at x = $(gridx_b[j_middle_column])")
        n = size(RX1, 1)
        for i in 1:n
            println("i = ", i, " y = ", gridy_vx[i], " RX1 = ", RX1[i, j_middle_column])
        end
        println()
        j_middle_column = size(RY1, 2) ÷ 2 + 1
        println("==> RY1 along column $j_middle_column at x = $(gridx_vy[j_middle_column])")
        n = size(RY1, 1)
        for i in 1:n
            println("i = ", i, " y = ", gridy_b[i], " RY1 = ", RY1[i, j_middle_column])
        end
        println()
    end
    return nothing
end

function print_viscoplastic_viscosity_etan0_etas0_along_middle(
    model::ModelData, 
    msg::String,
    print_row::Bool,
    print_column::Bool
)::Nothing
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array
    gridx_pr = model.grids.arrays.pressure.gridx_pr.array
    gridy_pr = model.grids.arrays.pressure.gridy_pr.array
    etan0 = model.stokes_continuity.arrays.viscosity.etan0.array
    etas0 = model.stokes_continuity.arrays.viscosity.etas0.array

    println("DEBUG: printing viscoplastic viscosity etas0 and etan0:", msg)

    ArrayUtils.print_min_max("etan0", etan0)
    ArrayUtils.print_min_max("etas0", etas0)

    if print_row
        i_middle_row = size(etan0, 1) ÷ 2 + 1
        println("==> Normal viscosity along row $i_middle_row at y = $(gridy_pr[i_middle_row])")
        n = size(etan0, 2)
        for j in 1:n
            println("j = ", j, " x = ", gridx_pr[j], " y = ", gridy_pr[i_middle_row], " etan = ", etan0[i_middle_row, j])
        end
        println()
        i_middle_row = size(etas0, 1) ÷ 2 + 1
        println("==> Shear viscosity along row $i_middle_row at y = $(gridy_pr[i_middle_row])")
        n = size(etas0, 2)
        for j in 1:n
            println("j = ", j, " x = ", gridx_b[j], " y = ", gridy_b[i_middle_row], " etas = ", etas0[i_middle_row, j])
        end
        println()
    end

    if print_column
        j_middle_column = size(etan0, 2) ÷ 2 + 1
        println("==> Normal viscosity along column $j_middle_column at x = $(gridx_pr[j_middle_column])")
        n = size(etan0, 1)
        for i in 1:n
            println("i = ", i, " y = ", gridy_pr[i], " etan0 = ", etan0[i, j_middle_column])
        end
        println()
        j_middle_column = size(etas0, 2) ÷ 2 + 1
        println("==> Shear viscosity along column $j_middle_column at x = $(gridx_b[j_middle_column])")
        n = size(etas0, 1)
        for i in 1:n
            println("i = ", i, " y = ", gridy_b[i], " etas0 = ", etas0[i, j_middle_column])
        end
        println()
    end
    return nothing
end

function print_stokes_system_of_equations_stats(
    Li::Vector{Int64},
    Lj::Vector{Int64},
    Lv::Vector{Float64},
    model::ModelData,
    msg::String
)::Nothing
    println("DEBUG: printing system of equations stats: ", msg)
    ArrayUtils.print_min_max("RHS", model.stokes_continuity.arrays.rhs.RHS.array)
    ArrayUtils.print_min_max("Li", Li)
    ArrayUtils.print_min_max("Lj", Lj)
    ArrayUtils.print_min_max("Lv", Lv)
    nnz = length(Lv)
    println("nnz = ", nnz)
    write_to_file = false
    if write_to_file
        ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
        # Write non-zero matrix values to file
        ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
        file_name_base = "large_mat_nz_jl_$(ntimestep)"
        fid_unique = check_file_name(file_name_base, 0)
        file_name = file_name_base * "_$(fid_unique).txt"
        println("=> file_name = ", file_name)
        open(file_name, "w") do f
            for i in 1:length(Lv)
                println(f, "$i\t", @sprintf("%.6e", Lv[i]))
            end
        end
        # Write RHS to file
        RHS = model.stokes_continuity.arrays.rhs.RHS.array
        file_name_base = "rhs_jl_$(ntimestep)"
        fid_unique = check_file_name(file_name_base, 0)
        file_name = file_name_base * "_$(fid_unique).txt"
        println("=> file_name = ", file_name)
        avg = 0.0
        nval = 0
        open(file_name, "w") do f
            n = length(RHS)
            for i in 1:n
                println(f, "$i\t", @sprintf("%.3e", RHS[i]))
                avg += RHS[i]
                nval += 1
            end
        end
        println("=> RHS avg = ", avg / nval)
    end
    return nothing
end

function print_heat_system_of_equations_stats(
    Li::Vector{Int64},
    Lj::Vector{Int64},
    Lv::Vector{Float64},
    model::ModelData,
    msg::String
)::Nothing
    println("DEBUG: printing system of equations stats for heat equation: ", msg)
    println("DEBUG: timestep_heat = ", model.timestep.parameters.thermal_loop.timestep_heat.value)
    ArrayUtils.print_min_max("tk0", model.heat_equation.arrays.temperature.tk0.array)
    ArrayUtils.print_min_max("rhocp1", model.heat_equation.arrays.rhocp.rhocp1.array)
    ArrayUtils.print_min_max("RT1", model.heat_equation.arrays.rhs.RT1.array)
    ArrayUtils.print_min_max("RHSheat", model.heat_equation.arrays.rhs.RHSheat.array)
    ArrayUtils.print_min_max("Li", Li)
    ArrayUtils.print_min_max("Lj", Lj)
    ArrayUtils.print_min_max("Lv", Lv)
    nnz = length(Lv)
    println("nnz = ", nnz)
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    # Write non-zero matrix values to file
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    file_name_base = "large_mat_heat_nz_jl_$(ntimestep)"
    fid_unique = check_file_name(file_name_base, 0)
    file_name = file_name_base * "_$(fid_unique).txt"
    println("=> file_name = ", file_name)
    open(file_name, "w") do f
        for i in 1:length(Lv)
            println(f, "$i\t", @sprintf("%.6e", Lv[i]))
        end
    end
    # Write RHS to file
    RHS = model.heat_equation.arrays.rhs.RHSheat.array
    file_name_base = "rhs_heat_jl_$(ntimestep)"
    fid_unique = check_file_name(file_name_base, 0)
    file_name = file_name_base * "_$(fid_unique).txt"
    println("=> file_name = ", file_name)
    avg = 0.0
    nval = 0
    open(file_name, "w") do f
        n = length(RHS)
        for i in 1:n
            println(f, "$i\t", @sprintf("%.3e", RHS[i]))
            avg += RHS[i]
            nval += 1
        end
    end
    println("=> RHSheat avg = ", avg / nval)
    return nothing
end

function check_file_name(file_name_base::String, fid::Int64)::Int64
    file_name = file_name_base * "_$(fid).txt"
    if isfile(file_name)
        fid = fid + 1
        fid = check_file_name(file_name_base, fid)
    end
    return fid
end

function print_stokes_solution_vector_stats(
    S::Vector{Float64},
    msg::String
)::Nothing
    println("DEBUG: printing stokes solution stats: ", msg)
    ArrayUtils.print_min_max("S", S)
    return nothing
end

function print_heat_solution_vector_stats(
    S::Vector{Float64},
    msg::String
)::Nothing
    println("DEBUG: printing heat solution stats: ", msg)
    ArrayUtils.print_min_max("S", S)
    return nothing
end

function print_stokes_solution_stats(model::ModelData, msg::String)::Nothing
    println("DEBUG: printing stokes solutions: ", msg)
    ArrayUtils.print_min_max("vx1", model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array)
    ArrayUtils.print_min_max("vy1", model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array)
    ArrayUtils.print_min_max("pr1", model.stokes_continuity.arrays.pressure.pr1.array)
    return nothing
end

function print_grid_strain_rate_and_spin_stats(model::ModelData, msg::String)::Nothing
    println("DEBUG: printing grid strain rate and spin stats: ", msg)
    ArrayUtils.print_min_max("exx", model.stokes_continuity.arrays.strain_rate_and_spin.exx.array)
    ArrayUtils.print_min_max("exy", model.stokes_continuity.arrays.strain_rate_and_spin.exy.array)
    ArrayUtils.print_min_max("esp", model.stokes_continuity.arrays.strain_rate_and_spin.esp.array)
    return nothing
end

function print_marker_stress_stats(model::ModelData, msg::String)::Nothing
    println("DEBUG: printing marker stress stats: ", msg)
    ArrayUtils.print_min_max("marker_sxy", model.markers.arrays.stress.marker_sxy.array)
    ArrayUtils.print_min_max("marker_sxx", model.markers.arrays.stress.marker_sxx.array)
    return nothing
end

function print_grid_stress_stats(model::ModelData, msg::String)::Nothing
    println("DEBUG: printing grid stress stats: ", msg)
    ArrayUtils.print_min_max("sxx2", model.stokes_continuity.arrays.stress.sxx2.array)
    ArrayUtils.print_min_max("sxy2", model.stokes_continuity.arrays.stress.sxy2.array)
    return nothing
end

function print_rhs_x_stokes_input_stats(model::ModelData, msg::String)::Nothing
    println("DEBUG: printing rhs_x_stokes input stats: ", msg)
    println("gravity_x = ", model.gravity.parameters.gravity_x.value)
    ArrayUtils.print_min_max("rho1", model.stokes_continuity.arrays.density.rho1.array)
    ArrayUtils.print_min_max("sxx0", model.stokes_continuity.arrays.stress.sxx0.array)
    ArrayUtils.print_min_max("sxy0", model.stokes_continuity.arrays.stress.sxy0.array)
    return nothing
end

function print_update_nodes_for_plastic_yielding_stats(
    etavp_new::Matrix{Float64},
    plastic_strain_rate::Matrix{Float64},
    plastic_yield_new::Matrix{Float64},
    nnode_yield::Int64,
    global_yield_error::Float64,
)::Nothing
    println("==> DEBUG: printing update nodes for plastic yielding stats")
    ArrayUtils.print_min_max("    > etavp_new", etavp_new)
    ArrayUtils.print_min_max("    > plastic_strain_rate", plastic_strain_rate)
    ArrayUtils.print_min_max("    > plastic_yield_new", plastic_yield_new)
    println("    > nnode_yield = ", nnode_yield)
    println("    > global_yield_error = ", global_yield_error)
    return nothing
end

function print_plasticity_loop_update_stats(model::ModelData)::Nothing
    etas1 = model.stokes_continuity.arrays.viscosity.etas1.array
    plastic_yield = model.stokes_continuity.arrays.plastic_def.plastic_yield.array
    plastic_strain_rate = model.stokes_continuity.arrays.strain_rate_and_spin.eii_plastic_basic.array
    dvxy_rel_L2 = model.stokes_continuity.parameters.solution_norms.dvxy_rel_L2.value
    println("==> DEBUG: printing plasticity loop update stats")
    ArrayUtils.print_min_max("    > etas1", etas1)
    ArrayUtils.print_min_max("    > plastic_yield", plastic_yield)
    ArrayUtils.print_min_max("    > plastic_strain_rate", plastic_strain_rate)
    println("    > dvxy_rel_L2 = ", dvxy_rel_L2)
    return nothing
end

function print_rhs_stats(model::ModelData, msg::String)::Nothing
    println("==> DEBUG: printing rhs stats: ", msg)
    ArrayUtils.print_min_max("    > RC1", model.stokes_continuity.arrays.rhs.RC1.array)
    ArrayUtils.print_min_max("    > RX1", model.stokes_continuity.arrays.rhs.RX1.array)
    ArrayUtils.print_min_max("    > RY1", model.stokes_continuity.arrays.rhs.RY1.array)
    ArrayUtils.print_min_max("    > RHS", model.stokes_continuity.arrays.rhs.RHS.array)
    ArrayUtils.print_min_max("    > rho1", model.stokes_continuity.arrays.density.rho1.array)
    ArrayUtils.print_min_max("    > sxx0", model.stokes_continuity.arrays.stress.sxx0.array)
    ArrayUtils.print_min_max("    > sxy0", model.stokes_continuity.arrays.stress.sxy0.array)
    return nothing
end

function print_sxy1_stats(model::ModelData, msg::String)::Nothing
    println("==> DEBUG: printing sxy1 stats: ", msg)
    ArrayUtils.print_min_max("    > sxy1", model.stokes_continuity.arrays.stress.sxy1.array)
    return nothing
end

function print_marker_sxy_stats(model::ModelData, msg::String)::Nothing
    println("==> DEBUG: printing marker sxy stats: ", msg)
    ArrayUtils.print_min_max("    > marker_sxy", model.markers.arrays.stress.marker_sxy.array)
    return nothing
end

end # module