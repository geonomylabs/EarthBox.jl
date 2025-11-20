module ArrayStats

using Printf
import Statistics
import ..LevelManager: LevelData, LevelData2d

function write_array_to_file(array::Array{Float64,3}, basename::String)::Nothing
    filename = "$basename"*"_jl.txt"
    open(filename, "w") do io
    for i in 1:size(array,1)
        for j in 1:size(array,2)
            for k in 1:size(array,3)
                    value = array[i,j,k]
                    println(io, "$i $j $k $(Printf.@sprintf("%.5e", value))")
                end
            end
        end
    end
    return nothing
end

function print_pre_vcycle_array_stats(level_vector::Vector{LevelData})
    print_density_stats(level_vector[1])
    print_viscosity_stats(level_vector[1])
    print_pressure_stats(level_vector[1])
    print_right_parts_stats(level_vector[1], "")
    levelnum = length(level_vector)
    for n = 1:levelnum
        println("")
        println(">> Viscosity interpolated to coarser levels")
        print_array_stats3d(level_vector[n].etan.array, "etan (level $(n))")
        print_array_stats3d(level_vector[n].etaxy.array, "etaxy (level $(n))")
        print_array_stats3d(level_vector[n].etaxz.array, "etaxz (level $(n))")
        print_array_stats3d(level_vector[n].etayz.array, "etayz (level $(n))")
        println(">> Backed up viscosity")
        print_array_stats3d(level_vector[n].etano.array, "etano (level $(n))")
        print_array_stats3d(level_vector[n].etaxyo.array, "etaxyo (level $(n))")
        print_array_stats3d(level_vector[n].etaxzo.array, "etaxzo (level $(n))")
        print_array_stats3d(level_vector[n].etayzo.array, "etayzo (level $(n))")
        println(">> Renormalized viscosity")
        print_array_stats3d(level_vector[n].etan.array, "etan (level $(n))")
        print_array_stats3d(level_vector[n].etaxy.array, "etaxy (level $(n))")
        print_array_stats3d(level_vector[n].etaxz.array, "etaxz (level $(n))")
        print_array_stats3d(level_vector[n].etayz.array, "etayz (level $(n))")
    end
end

function print_array_stats3d(array::Union{Array{Float64, 3}, Array{Int64, 3}, Vector{Float64}}, msg::String)
    # Print min, max and average values
    min_value = minimum(array)
    max_value = maximum(array)
    avg_value = Statistics.mean(array)
    println("$msg min: $min_value, max: $max_value, avg: $avg_value")
end

function print_array_stats(
    array::Union{Array{Float64, 3}, Array{Int64, 3}, Matrix{Float64}, Matrix{Int64}, Vector{Float64}}, 
    msg::String
)::Nothing
    # Print min, max and average values
    min_value = minimum(array)
    max_value = maximum(array)
    avg_value = Statistics.mean(array)
    println("$msg min: $min_value, max: $max_value, avg: $avg_value")
    return nothing
end

function print_density_stats(level_data::Union{LevelData, LevelData2d})
    # Print min, max and average density values
    min_density = minimum(level_data.rho.array)
    max_density = maximum(level_data.rho.array)
    avg_density = Statistics.mean(level_data.rho.array)
    println("rho: min: $min_density, max: $max_density, avg: $avg_density")
end

function print_viscosity_stats(level_data::LevelData)
    # Print min, max and average viscosity values
    min_viscosity = minimum(level_data.etan.array)
    max_viscosity = maximum(level_data.etan.array)
    avg_viscosity = Statistics.mean(level_data.etan.array)
    println("etan: min: $min_viscosity, max: $max_viscosity, avg: $avg_viscosity")

    min_viscosity = minimum(level_data.etaxy.array)
    max_viscosity = maximum(level_data.etaxy.array)
    avg_viscosity = Statistics.mean(level_data.etaxy.array)
    println("etaxy: min: $min_viscosity, max: $max_viscosity, avg: $avg_viscosity")
    
    min_viscosity = minimum(level_data.etaxz.array)
    max_viscosity = maximum(level_data.etaxz.array)
    avg_viscosity = Statistics.mean(level_data.etaxz.array)
    println("etaxz: min: $min_viscosity, max: $max_viscosity, avg: $avg_viscosity")
    
    min_viscosity = minimum(level_data.etayz.array)
    max_viscosity = maximum(level_data.etayz.array)
    avg_viscosity = Statistics.mean(level_data.etayz.array)
    println("etayz: min: $min_viscosity, max: $max_viscosity, avg: $avg_viscosity")
end

function print_viscosity_stats(level_data::LevelData2d)
    # Print min, max and average viscosity values
    min_viscosity = minimum(level_data.etan.array)
    max_viscosity = maximum(level_data.etan.array)
    avg_viscosity = Statistics.mean(level_data.etan.array)
    println("etan: min: $min_viscosity, max: $max_viscosity, avg: $avg_viscosity")

    min_viscosity = minimum(level_data.etas.array)
    max_viscosity = maximum(level_data.etas.array)
    avg_viscosity = Statistics.mean(level_data.etas.array)
    println("etas: min: $min_viscosity, max: $max_viscosity, avg: $avg_viscosity")
end

function print_pressure_stats(level_data::Union{LevelData, LevelData2d})
    # Print min, max and average pressure values
    min_pressure = minimum(level_data.pr.array)
    max_pressure = maximum(level_data.pr.array)
    avg_pressure = Statistics.mean(level_data.pr.array)
    println("pr: min: $min_pressure, max: $max_pressure, avg: $avg_pressure")
end

function print_solution_stats(
    level_data::LevelData,
    msg::String
)::Nothing
    println(">>>> $msg solution statistics:")
    min_vx = minimum(level_data.vx.array)
    max_vx = maximum(level_data.vx.array)
    avg_vx = Statistics.mean(level_data.vx.array)
    min_vy = minimum(level_data.vy.array)
    max_vy = maximum(level_data.vy.array)
    avg_vy = Statistics.mean(level_data.vy.array)
    min_vz = minimum(level_data.vz.array)
    max_vz = maximum(level_data.vz.array)
    avg_vz = Statistics.mean(level_data.vz.array)
    min_pr = minimum(level_data.pr.array)
    max_pr = maximum(level_data.pr.array)
    avg_pr = Statistics.mean(level_data.pr.array)
    println("vx_min: $min_vx, vx_max: $max_vx, vx_avg: $avg_vx")
    println("vy_min: $min_vy, vy_max: $max_vy, vy_avg: $avg_vy")
    println("vz_min: $min_vz, vz_max: $max_vz, vz_avg: $avg_vz")
    println("pr_min: $min_pr, pr_max: $max_pr, pr_avg: $avg_pr")
end

function print_solution_stats(
    level_data::LevelData2d,
    msg::String
)::Nothing
    println(">>>> $msg solution statistics:")
    min_vx = minimum(level_data.vx.array)
    max_vx = maximum(level_data.vx.array)
    avg_vx = Statistics.mean(level_data.vx.array)
    min_vy = minimum(level_data.vy.array)
    max_vy = maximum(level_data.vy.array)
    avg_vy = Statistics.mean(level_data.vy.array)
    min_pr = minimum(level_data.pr.array)
    max_pr = maximum(level_data.pr.array)
    avg_pr = Statistics.mean(level_data.pr.array)
    println("vx_min: $min_vx, vx_max: $max_vx, vx_avg: $avg_vx")
    println("vy_min: $min_vy, vy_max: $max_vy, vy_avg: $avg_vy")
    println("pr_min: $min_pr, pr_max: $max_pr, pr_avg: $avg_pr")
end

function print_residual_stats(
    resx::Array{Float64,3},
    resy::Array{Float64,3},
    resz::Array{Float64,3},
    resc::Array{Float64,3},
    msg::String
)::Nothing
    println(">>>> $msg residual statistics:")
    min_resx = minimum(resx)
    max_resx = maximum(resx)
    avg_resx = Statistics.mean(resx)
    min_resy = minimum(resy)
    max_resy = maximum(resy)
    avg_resy = Statistics.mean(resy)
    min_resz = minimum(resz)
    max_resz = maximum(resz)
    avg_resz = Statistics.mean(resz)
    min_resc = minimum(resc)
    max_resc = maximum(resc)
    avg_resc = Statistics.mean(resc)
    println("resx_min: $min_resx, resx_max: $max_resx, resx_avg: $avg_resx")
    println("resy_min: $min_resy, resy_max: $max_resy, resy_avg: $avg_resy")
    println("resz_min: $min_resz, resz_max: $max_resz, resz_avg: $avg_resz")
    println("resc_min: $min_resc, resc_max: $max_resc, resc_avg: $avg_resc")
end

function print_residual_stats_2d(
    resx::Array{Float64,2},
    resy::Array{Float64,2},
    resc::Array{Float64,2},
    msg::String
)::Nothing
    println(">>>> $msg residual statistics:")
    min_resx = minimum(resx)
    max_resx = maximum(resx)
    avg_resx = Statistics.mean(resx)
    min_resy = minimum(resy)
    max_resy = maximum(resy)
    avg_resy = Statistics.mean(resy)
    min_resc = minimum(resc)
    max_resc = maximum(resc)
    avg_resc = Statistics.mean(resc)
    println("resx_min: $min_resx, resx_max: $max_resx, resx_avg: $avg_resx")
    println("resy_min: $min_resy, resy_max: $max_resy, resy_avg: $avg_resy")
    println("resc_min: $min_resc, resc_max: $max_resc, resc_avg: $avg_resc")
end

function print_right_parts_stats(
    level_data::LevelData,
    msg::String
)::Nothing
    println(">>>> $msg right parts statistics:")
    min_rx = minimum(level_data.RX.array)
    max_rx = maximum(level_data.RX.array)
    avg_rx = Statistics.mean(level_data.RX.array)
    min_ry = minimum(level_data.RY.array)
    max_ry = maximum(level_data.RY.array)
    avg_ry = Statistics.mean(level_data.RY.array)
    min_rz = minimum(level_data.RZ.array)
    max_rz = maximum(level_data.RZ.array)
    avg_rz = Statistics.mean(level_data.RZ.array)
    min_rc = minimum(level_data.RC.array)
    max_rc = maximum(level_data.RC.array)
    avg_rc = Statistics.mean(level_data.RC.array)
    println("rx_min: $min_rx, rx_max: $max_rx, rx_avg: $avg_rx")
    println("ry_min: $min_ry, ry_max: $max_ry, ry_avg: $avg_ry")
    println("rz_min: $min_rz, rz_max: $max_rz, rz_avg: $avg_rz")
    println("rc_min: $min_rc, rc_max: $max_rc, rc_avg: $avg_rc")
end

function print_right_parts_stats(
    level_data::LevelData2d,
    msg::String
)::Nothing
    println(">>>> $msg right parts statistics:")
    min_rx = minimum(level_data.RX.array)
    max_rx = maximum(level_data.RX.array)
    avg_rx = Statistics.mean(level_data.RX.array)
    min_ry = minimum(level_data.RY.array)
    max_ry = maximum(level_data.RY.array)
    avg_ry = Statistics.mean(level_data.RY.array)
    min_rc = minimum(level_data.RC.array)
    max_rc = maximum(level_data.RC.array)
    avg_rc = Statistics.mean(level_data.RC.array)
    println("rx_min: $min_rx, rx_max: $max_rx, rx_avg: $avg_rx")
    println("ry_min: $min_ry, ry_max: $max_ry, ry_avg: $avg_ry")
    println("rc_min: $min_rc, rc_max: $max_rc, rc_avg: $avg_rc")
end

function print_correction_stats(
    dvx::Array{Float64,3},
    dvy::Array{Float64,3},
    dvz::Array{Float64,3},
    dpr::Array{Float64,3},
    msg::String
)::Nothing
    println(">>>> $msg correction statistics:")
    min_dvx = minimum(dvx)
    max_dvx = maximum(dvx)
    avg_dvx = Statistics.mean(dvx)
    min_dvy = minimum(dvy)
    max_dvy = maximum(dvy)
    avg_dvy = Statistics.mean(dvy)
    min_dvz = minimum(dvz)
    max_dvz = maximum(dvz)
    avg_dvz = Statistics.mean(dvz)
    min_dpr = minimum(dpr)
    max_dpr = maximum(dpr)
    avg_dpr = Statistics.mean(dpr)
    println("dvx_min: $min_dvx, dvx_max: $max_dvx, dvx_avg: $avg_dvx")
    println("dvy_min: $min_dvy, dvy_max: $max_dvy, dvy_avg: $avg_dvy")
    println("dvz_min: $min_dvz, dvz_max: $max_dvz, dvz_avg: $avg_dvz")
    println("dpr_min: $min_dpr, dpr_max: $max_dpr, dpr_avg: $avg_dpr")
end

function print_correction_stats_2d(
    dvx::Array{Float64,2},
    dvy::Array{Float64,2},
    dpr::Array{Float64,2},
    msg::String
)::Nothing
    println(">>>> $msg correction statistics:")
    min_dvx = minimum(dvx)
    max_dvx = maximum(dvx)
    avg_dvx = Statistics.mean(dvx)
    min_dvy = minimum(dvy)
    max_dvy = maximum(dvy)
    avg_dvy = Statistics.mean(dvy)
    min_dpr = minimum(dpr)
    max_dpr = maximum(dpr)
    avg_dpr = Statistics.mean(dpr)
    println("dvx_min: $min_dvx, dvx_max: $max_dvx, dvx_avg: $avg_dvx")
    println("dvy_min: $min_dvy, dvy_max: $max_dvy, dvy_avg: $avg_dvy")
    println("dpr_min: $min_dpr, dpr_max: $max_dpr, dpr_avg: $avg_dpr")
end

end # module