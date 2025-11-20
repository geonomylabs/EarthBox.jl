"""
    Smoother

This module provides a function that makes specified number (smoothing_iterations) of Gauss Seidel 
iterations using relaxation coefficients (krelaxs, krelaxc) for Stokes and Continuity 
equations defined on 3D staggered grid with specified resolution (xnum, ynum, znum) 
and grid steps (xstp, ystp, zstp).

Given distribution of right parts for all equations (RX, RY, RZ, RC) on the grid 
and given variable viscosity (etaxy, etaxz, etayz, etan) on the grid.
Pressure is normalized relative to given value (prnorm) in the first cell.

The function returns new approximation for velocity and pressure (vx, vy, vz, pr)
and distribution of residuals (resx, resy, resz, resc).

# Staggered Grid for Multigrid
# 
#     vx       vx       vx    
#
# vy  +---vy---+---vy---+   vy
#     |        |        |
#     vx   P   vx   P   vx    
#     |        |        |
# vy  +---vy---+---vy---+   vy
#     |        |        |
#     vx   P   vx   P   vx    
#     |        |        |
# vy  +---vy---+---vy---+   vy
#
#     vx       vx       vx    
# 
# Lines show basic grid
# Basic (density) nodes are shown with +
# Ghost nodes shown outside the basic grid
# are used for boundary conditions
"""
module Smoother

include("residuals/Residuals.jl")
include("solver/StokesSolverGaussSeidel.jl")
include("pressure_update/PressureUpdate.jl")

import ..MultigridDataManager: MultigridData3d, MultigridData2d
import ..MultigridDataManager.MultigridStructures: RelaxationParameters
import ..BoundaryConditions: set_boundary_conditions3d!
import ..BoundaryConditions: set_boundary_conditions2d!
import ..LevelManager: LevelData, LevelData2d
import .StokesSolverGaussSeidel: solve_stokes_equations_gauss_seidel!
import .PressureUpdate: update_pressure!
import .Residuals: compute_residuals!

"""
    stokes_continuity3d_viscous_smoother!(pressure_bc, level_data, smoothing_iterations, relaxation)

Perform Gauss-Seidel iterations for Stokes and Continuity equations on a 3D staggered grid.

# Arguments
- `pressure_bc`: Pressure boundary condition
- `level_data`: Level data
- `smoothing_iterations`: Number of Gauss-Seidel smoothing iterations
- `relaxation`: Relaxation parameters

# Returns
- `resx`, `resy`, `resz`, `resc`: Residuals for all equations
"""
function stokes_continuity3d_viscous_smoother!(
    pressure_bc::Float64,
    level_data::LevelData,
    smoothing_iterations::Vector{Int64},
    relaxation::RelaxationParameters
)::Tuple{Array{Float64,3}, Array{Float64,3}, Array{Float64,3}, Array{Float64,3}}
    solve_stokes_equations_gauss_seidel!(level_data, smoothing_iterations, relaxation)
    update_pressure!(pressure_bc, level_data)
    resx, resy, resz, resc = compute_residuals!(level_data)
    return resx, resy, resz, resc
end

function stokes_continuity2d_viscous_smoother!(
    pressure_bc::Float64,
    level_data::LevelData2d,
    smoothing_iterations::Vector{Int64},
    relaxation::RelaxationParameters
)::Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}}
    solve_stokes_equations_gauss_seidel!(level_data, smoothing_iterations, relaxation)
    update_pressure!(pressure_bc, level_data)
    resx, resy, resc = compute_residuals!(level_data)
    return resx, resy, resc
end

function update_rhs_parts_on_level1_using_global_level0_residuals(
    multigrid_data::MultigridData3d
)::Tuple{Array{Float64,3}, Array{Float64,3}, Array{Float64,3}, Array{Float64,3}}
    level0 = multigrid_data.level0
    level1 = multigrid_data.level_vector[1]
    pressure_bc = multigrid_data.pressure_bc
    (
        ΔRx_L0, ΔRy_L0, ΔRz_L0, ΔRc_L0
    ) = calculate_residuals_without_smoothing(pressure_bc, level0)
    level1.RX.array .= ΔRx_L0
    level1.RY.array .= ΔRy_L0
    level1.RZ.array .= ΔRz_L0
    level1.RC.array .= ΔRc_L0
    return ΔRx_L0, ΔRy_L0, ΔRz_L0, ΔRc_L0
end

function update_rhs_parts_on_level1_using_global_level0_residuals(
    multigrid_data::MultigridData2d
)::Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}}
    level0 = multigrid_data.level0
    level1 = multigrid_data.level_vector[1]
    pressure_bc = multigrid_data.pressure_bc
    (
        ΔRx_L0, ΔRy_L0, ΔRc_L0
    ) = calculate_residuals_without_smoothing(pressure_bc, level0)
    level1.RX.array .= ΔRx_L0
    level1.RY.array .= ΔRy_L0
    level1.RC.array .= ΔRc_L0
    return ΔRx_L0, ΔRy_L0, ΔRc_L0
end

function calculate_residuals_without_smoothing(
    pressure_bc::Float64,
    level_data::LevelData,
)::Tuple{Array{Float64,3}, Array{Float64,3}, Array{Float64,3}, Array{Float64,3}}
    set_boundary_conditions3d!(level_data)
    update_pressure!(pressure_bc, level_data)
    resx, resy, resz, resc = compute_residuals!(level_data)
    return resx, resy, resz, resc
end

function calculate_residuals_without_smoothing(
    pressure_bc::Float64,
    level_data::LevelData2d,
)::Tuple{Array{Float64,2}, Array{Float64,2}, Array{Float64,2}}
    set_boundary_conditions2d!(level_data)
    update_pressure!(pressure_bc, level_data)
    resx, resy, resc = compute_residuals!(level_data)
    return resx, resy, resc
end

end # module 