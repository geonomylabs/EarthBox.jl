module SolveStokes3dOpt

import ....MultigridDataManager.MultigridStructures: RelaxationParameters
import ....LevelManager: LevelData
import ....Domain: on_vx_boundary3d, on_vy_boundary3d, on_vz_boundary3d
import ....ArrayStats
import ...Residuals

function solve_stokes_continuity_equations3d!(
    relaxation::RelaxationParameters,
    level_data::LevelData
)::Nothing
    xnum = level_data.grid.parameters.geometry.xnum.value
    ynum = level_data.grid.parameters.geometry.ynum.value
    znum = level_data.grid.parameters.geometry.znum.value

    Θ_stokes = relaxation.relax_stokes
    Θ_continuity = relaxation.relax_continuity

    Threads.@threads for k = 1:znum+1
        for j = 1:xnum+1
            for i = 1:ynum+1
                if j < xnum+1 # vx is not defined beyond xnum
                    if !on_vx_boundary3d(i, j, k, ynum, xnum, znum)
                        update_vx!(i, j, k, Θ_stokes, level_data)
                    end
                end
                if i < ynum+1 # vy is not defined beyond ynum
                    if !on_vy_boundary3d(i, j, k, ynum, xnum, znum)
                        update_vy!(i, j, k, Θ_stokes, level_data)
                    end
                end
                if k < znum+1 # vz is not defined beyond znum
                    if !on_vz_boundary3d(i, j, k, ynum, xnum, znum)
                        update_vz!(i, j, k, Θ_stokes, level_data)
                    end
                end
                if i < ynum && j < xnum && k < znum
                    update_pressure!(i, j, k, Θ_continuity, level_data)
                end
            end
        end
    end
    return nothing
end

# x-Stokes equation dSIGMAxx/dx+dSIGMAxy/dy+SIGMAxz/dz-dP/dx=RX
function update_vx!(
    i::Int64,
    j::Int64,
    k::Int64,
    Θ_stokes::Float64,
    level_data::LevelData
)::Nothing
    ΔR, Coef_vxC = Residuals.calculate_vx_residual(i, j, k, level_data)
    level_data.vx.array[i,j,k] += ΔR/Coef_vxC*Θ_stokes
    return nothing
end

# y-Stokes equation dSIGMAyx/dx+dSIGMAyy/dy+SIGMAyz/dz-dP/dy=RY
function update_vy!(
    i::Int64,
    j::Int64,
    k::Int64,
    Θ_stokes::Float64,
    level_data::LevelData
)::Nothing
    ΔR, Coef_vyC = Residuals.calculate_vy_residual(i, j, k, level_data)
    level_data.vy.array[i,j,k] += ΔR/Coef_vyC*Θ_stokes
    return nothing
end

# z-Stokes equation dSIGMAzx/dx+dSIGMAzy/dy+SIGMAzz/dz-dP/dz=RZ
function update_vz!(
    i::Int64,
    j::Int64,
    k::Int64,
    Θ_stokes::Float64,
    level_data::LevelData
)::Nothing
    ΔR, Coef_vzC = Residuals.calculate_vz_residual(i, j, k, level_data)
    level_data.vz.array[i,j,k] += ΔR/Coef_vzC*Θ_stokes
    return nothing
end

function update_pressure!(
    i::Int64,
    j::Int64,
    k::Int64,
    Θ_continuity::Float64,
    level_data::LevelData
)::Nothing
    ΔR = Residuals.calculate_pressure_residual(i, j, k, level_data)
    level_data.pr.array[i,j,k] += level_data.etan.array[i,j,k]*ΔR*Θ_continuity
    return nothing
end

end # module