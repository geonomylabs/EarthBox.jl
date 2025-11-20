module SmoothAndProlongate

import ...MultigridDataManager: MultigridData3d, MultigridData2d
import ...Smoother: stokes_continuity3d_viscous_smoother!
import ...Smoother: stokes_continuity2d_viscous_smoother!
import ...Prolongation: prolongate_stokes3d_solution
import ...Prolongation: prolongate_stokes2d_solution

function smooth_and_prolongate!(
    multigrid_data::MultigridData3d,
)::Tuple{Array{Float64,3}, Array{Float64,3}, Array{Float64,3}, Array{Float64,3}}
    pressure_bc = multigrid_data.pressure_bc
    smoothing_iterations = multigrid_data.smoothing_iterations
    relaxation = multigrid_data.relaxation
    level_vector = multigrid_data.level_vector

    levelnum = length(level_vector)
    ΔRxL1 = zeros(Float64, 2, 2, 2)
    ΔRyL1 = zeros(Float64, 2, 2, 2)
    ΔRzL1 = zeros(Float64, 2, 2, 2)
    ΔRcL1 = zeros(Float64, 2, 2, 2)
    for n = levelnum:-1:1
        if n == 1
            ΔRxL1, ΔRyL1, ΔRzL1, ΔRcL1 = stokes_continuity3d_viscous_smoother!(
                    pressure_bc, level_vector[n], smoothing_iterations, relaxation)
        else
            _, _, _, _ = stokes_continuity3d_viscous_smoother!(
                0.0, level_vector[n], smoothing_iterations, relaxation)
            dvx, dvy, dvz, dpr = prolongate_stokes3d_solution(n, level_vector)
            relax_velocity = multigrid_data.relaxation.relax_velocity
            relax_pressure = multigrid_data.relaxation.relax_pressure
            level_vector[n-1].vx.array .+= dvx * relax_velocity
            level_vector[n-1].vy.array .+= dvy * relax_velocity
            level_vector[n-1].vz.array .+= dvz * relax_velocity
            level_vector[n-1].pr.array .+= dpr * relax_pressure
        end
    end
    return ΔRxL1, ΔRyL1, ΔRzL1, ΔRcL1
end

function smooth_and_prolongate!(
    multigrid_data::MultigridData2d,
)::Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}}
    pressure_bc = multigrid_data.pressure_bc
    smoothing_iterations = multigrid_data.smoothing_iterations
    relaxation = multigrid_data.relaxation
    level_vector = multigrid_data.level_vector

    levelnum = length(level_vector)
    ΔRxL1 = zeros(Float64, 2, 2)
    ΔRyL1 = zeros(Float64, 2, 2)
    ΔRcL1 = zeros(Float64, 2, 2)
    for n = levelnum:-1:1
        if n == 1
            ΔRxL1, ΔRyL1, ΔRcL1 = stokes_continuity2d_viscous_smoother!(
                    pressure_bc, level_vector[n], smoothing_iterations, relaxation)
        else
            _, _, _ = stokes_continuity2d_viscous_smoother!(
                0.0, level_vector[n], smoothing_iterations, relaxation)
            dvx, dvy, dpr = prolongate_stokes2d_solution(n, level_vector)
            relax_velocity = multigrid_data.relaxation.relax_velocity
            relax_pressure = multigrid_data.relaxation.relax_pressure
            level_vector[n-1].vx.array .+= dvx * relax_velocity
            level_vector[n-1].vy.array .+= dvy * relax_velocity
            level_vector[n-1].pr.array .+= dpr * relax_pressure
        end
    end
    return ΔRxL1, ΔRyL1, ΔRcL1
end

end