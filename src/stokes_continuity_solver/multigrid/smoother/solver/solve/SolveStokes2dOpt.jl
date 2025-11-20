module SolveStokes2dOpt

import ....MultigridDataManager.MultigridStructures: RelaxationParameters
import ....LevelManager: LevelData2d
import ....Domain: on_vx_boundary2d, on_vy_boundary2d
import ....ArrayStats
import ...Residuals

function solve_stokes_continuity_equations2d!(
    relaxation::RelaxationParameters,
    level_data::LevelData2d
)::Nothing
    xnum = level_data.grid.parameters.geometry.xnum.value
    ynum = level_data.grid.parameters.geometry.ynum.value

    Θ_stokes = relaxation.relax_stokes
    Θ_continuity = relaxation.relax_continuity

    Threads.@threads for j = 1:xnum+1
        for i = 1:ynum+1
            if j < xnum+1 # vx is not defined beyond xnum
                if !on_vx_boundary2d(i, j, ynum, xnum)
                    update_vx!(i, j, Θ_stokes, level_data)
                end
            end
            if i < ynum+1 # vy is not defined beyond ynum
                if !on_vy_boundary2d(i, j, ynum, xnum)
                    update_vy!(i, j, Θ_stokes, level_data)
                end
            end
            if i < ynum && j < xnum
                update_pressure!(i, j, Θ_continuity, level_data)
            end
        end
    end
    return nothing
end

# x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy+SIGMAxz/dz-dP/dx=RX
function update_vx!(
    i::Int64,
    j::Int64,
    Θ_stokes::Float64,
    level_data::LevelData2d
)::Nothing
    ΔR, Coef_vxC = Residuals.calculate_vx_residual(i, j, level_data)
    level_data.vx.array[i,j] += ΔR/Coef_vxC*Θ_stokes
    return nothing
end

# y-Stokes equation dSIGMAyx/dx+dSIGMAyy/dy+SIGMAyz/dz-dP/dy=RY
function update_vy!(
    i::Int64,
    j::Int64,
    Θ_stokes::Float64,
    level_data::LevelData2d
)::Nothing
    ΔR, Coef_vyC = Residuals.calculate_vy_residual(i, j, level_data)
    level_data.vy.array[i,j] += ΔR/Coef_vyC*Θ_stokes
    return nothing
end

function update_pressure!(
    i::Int64,
    j::Int64,
    Θ_continuity::Float64,
    level_data::LevelData2d
)::Nothing
    ΔR = Residuals.calculate_pressure_residual(i, j, level_data)
    level_data.pr.array[i,j] += level_data.etan.array[i,j]*ΔR*Θ_continuity
    return nothing
end

end # module