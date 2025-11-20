"""
    ResidualPlotting

Module for plotting residual surfaces for Stokes and continuity equations in 3D.
This module provides functions to create surface plots similar to the MATLAB 
implementation in Variable_viscosity3D_Multigrid.m
"""
module ResidualPlotting

using Plots
import CairoMakie
import EarthBox.ModelDataContainer.MultiGrids2dContainer: MultigridData
import ..LevelManager: LevelData
import ..LevelManager: LevelData2d

const OUTPUT_DIR_NAME = "output_multigrid_test"

"""
    plot_residual_surfaces(
        resx0::Array{Float64, 3},
        resy0::Array{Float64, 3}, 
        resz0::Array{Float64, 3},
        resc0::Array{Float64, 3},
        znum::Int64,
        ivcycle::Int64 = 0
    )

Plot residual surfaces for x-Stokes, y-Stokes, z-Stokes, and continuity equations.
This function replicates the MATLAB surface plotting functionality from the original code.

# Arguments
- `resx0`: 3D array of x-Stokes residuals
- `resy0`: 3D array of y-Stokes residuals  
- `resz0`: 3D array of z-Stokes residuals
- `resc0`: 3D array of continuity residuals
- `znum`: Grid dimension
- `ivcycle`: Optional iteration number for title (default: 0)

# Returns
- A plot object containing the 4 subplots
"""
function plot_residual_surfaces(
    resx0::Array{Float64, 3},
    resy0::Array{Float64, 3}, 
    resz0::Array{Float64, 3},
    resc0::Array{Float64, 3},
    level_data::LevelData,
    ivcycle::Int64 = 0;
    msg::String="MainCycle"
)::Nothing
    min_delta = 1e-15
    # Calculate the middle z-slice index
    znum = level_data.grid.parameters.geometry.znum.value
    vy = level_data.vy.array

    z_mid = div(znum - 1, 2) + 1
    
    # Create coordinate grids for each residual array separately
    y_coords_x = 1:size(resx0, 1)
    x_coords_x = 1:size(resx0, 2)
    resx0_min = minimum(resx0)
    resx0_max = maximum(resx0)
    delta = resx0_max - resx0_min
    #println("resx0_min = $resx0_min, resx0_max = $resx0_max, delta = $delta")
    if delta < min_delta
        resx0_max = resx0_min + min_delta
    end
    p1 = surface(
        x_coords_x, y_coords_x, resx0[:, :, z_mid],
        title = ivcycle > 0 ? "x-Stokes residual, cycle = $ivcycle" : "x-Stokes residual",
        zlabel = "residual x-Stokes", 
        color = :viridis,
        camera = (45, 45),
        colorbar = false,
        zlims = (resx0_min, resx0_max)
        )
    
    y_coords_y = 1:size(resy0, 1)
    x_coords_y = 1:size(resy0, 2)
    resy0_min = minimum(resy0)
    resy0_max = maximum(resy0)
    delta = resy0_max - resy0_min
    #println("resy0_min = $resy0_min, resy0_max = $resy0_max, delta = $delta")
    if delta < min_delta
        resy0_max = resy0_min + min_delta
    end
    p2 = surface(
        x_coords_y, y_coords_y, resy0[:, :, z_mid],
        title = ivcycle > 0 ? "y-Stokes residual, cycle = $ivcycle" : "y-Stokes residual", 
        zlabel = "residual y-Stokes",
        color = :viridis,
        camera = (45, 45),
        colorbar = false,
        zlims = (resy0_min, resy0_max)
        )
    
    y_coords_z = 1:size(resz0, 1)
    x_coords_z = 1:size(resz0, 2)
    resz0_min = minimum(resz0)
    resz0_max = maximum(resz0)
    delta = resz0_max - resz0_min
    #println("resz0_min = $resz0_min, resz0_max = $resz0_max, delta = $delta")
    if delta < min_delta
        resz0_max = resz0_min + min_delta
    end
    p3 = surface(
        x_coords_z, y_coords_z, resz0[:, :, z_mid],
        title = ivcycle > 0 ? "z-Stokes residual, cycle = $ivcycle" : "z-Stokes residual",
        zlabel = "residual z-Stokes", 
        color = :viridis,
        camera = (45, 45),
        colorbar = false,
        zlims = (resz0_min, resz0_max)
        )
    
    y_coords_c = 1:size(resc0, 1)
    x_coords_c = 1:size(resc0, 2)
    resc0_min = minimum(resc0)
    resc0_max = maximum(resc0)
    delta = resc0_max - resc0_min
    #println("resc0_min = $resc0_min, resc0_max = $resc0_max, delta = $delta")
    if delta < min_delta
        resc0_max = resc0_min + min_delta
    end
    p4 = surface(
        x_coords_c, y_coords_c, resc0[:, :, z_mid],
        title = ivcycle > 0 ? "Continuity residual, cycle = $ivcycle" : "Continuity residual",
        zlabel = "residual continuity",
        color = :viridis,
        camera = (45, 45),
        colorbar = false,
        zlims = (resc0_min, resc0_max)
        )

    
    y_coords_c = 1:size(vy, 1)
    x_coords_c = 1:size(vy, 2)
    vy_min = minimum(vy)
    vy_max = maximum(vy)
    delta = vy_max - vy_min
    #println("vy_min = $vy_min, vy_max = $vy_max, delta = $delta")
    if delta < min_delta
        vy_max = vy_min + min_delta
    end
    p5 = surface(
        x_coords_c, y_coords_c, vy[:, :, z_mid],
        title = "Vy at i = $z_mid (z_mid)",
        zlabel = "Vy (m/s)",
        color = :viridis,
        camera = (45, 45),
        colorbar = false,
        zlims = (vy_min, vy_max)
        )
    
    # Combine all plots into a single figure and save
    p = plot(
        p1, p2, p3, p4, p5,
        layout = (2, 3),
        size = (1200, 800), # Increased size to give more room
        margin = 5Plots.mm, # Add margin between plots
        plot_title = "Residual Surfaces $msg"
        )
    # Check if output directory exists, if not create it
    check_dir(OUTPUT_DIR_NAME)
    savefig(p, joinpath(OUTPUT_DIR_NAME, "residual_surfaces_iter$(ivcycle).png"))
    cleanup()
    return nothing
end

function plot_residual_surfaces(
    resx0::Array{Float64, 2},
    resy0::Array{Float64, 2}, 
    resc0::Array{Float64, 2},
    level_data::LevelData2d,
    ivcycle::Int64=0;
    msg::String="MainCycle"
)::Nothing
    min_delta = 1e-15

    # Calculate the middle z-slice index
    vy = level_data.vy.array

    # Create coordinate grids for each residual array separately

    y_coords_x = 1:size(resx0, 1)
    x_coords_x = 1:size(resx0, 2)
    resx0_min = minimum(resx0)
    resx0_max = maximum(resx0)
    delta = resx0_max - resx0_min
    #println("resx0_min = $resx0_min, resx0_max = $resx0_max, delta = $delta")
    if delta < min_delta
        resx0_max = resx0_min + min_delta
    end
    p1 = surface(
        x_coords_x, y_coords_x, resx0,
        title = ivcycle > 0 ? "x-Stokes residual, cycle = $ivcycle" : "x-Stokes residual",
        zlabel = "residual x-Stokes", 
        color = :viridis,
        camera = (45, 45),
        colorbar = false,
        zlims = (resx0_min, resx0_max),
        reuse=false
        )
    y_coords_y = 1:size(resy0, 1)
    x_coords_y = 1:size(resy0, 2)
    resy0_min = minimum(resy0)
    resy0_max = maximum(resy0)
    delta = resy0_max - resy0_min
    #println("resy0_min = $resy0_min, resy0_max = $resy0_max, delta = $delta")
    if delta < min_delta
        resy0_max = resy0_min + min_delta
    end
    p2 = surface(
        x_coords_y, y_coords_y, resy0,
        title = ivcycle > 0 ? "y-Stokes residual, cycle = $ivcycle" : "y-Stokes residual", 
        zlabel = "residual y-Stokes",
        color = :viridis,
        camera = (45, 45),
        colorbar = false,
        zlims = (resy0_min, resy0_max),
        reuse=false
        )
    y_coords_c = 1:size(resc0, 1)
    x_coords_c = 1:size(resc0, 2)
    resc0_min = minimum(resc0)
    resc0_max = maximum(resc0)
    delta = resc0_max - resc0_min
    #println("resc0_min = $resc0_min, resc0_max = $resc0_max, delta = $delta")
    if delta < min_delta
        resc0_max = resc0_min + min_delta
    end
    p4 = surface(
        x_coords_c, y_coords_c, resc0,
        title = ivcycle > 0 ? "Continuity residual, cycle = $ivcycle" : "Continuity residual",
        zlabel = "residual continuity",
        color = :viridis,
        camera = (45, 45),
        colorbar = false,
        zlims = (resc0_min, resc0_max),
        reuse=false
        )

    y_coords_c = 1:size(vy, 1)
    x_coords_c = 1:size(vy, 2)
    vy_min = minimum(vy)
    vy_max = maximum(vy)
    delta = vy_max - vy_min
    #println("vy_min = $vy_min, vy_max = $vy_max, delta = $delta")
    if delta < min_delta
        vy_max = vy_min + min_delta
    end
    p5 = surface(
        x_coords_c, y_coords_c, vy,
        title = "Vy cycle = $ivcycle",
        zlabel = "Vy (m/s)", 
        color = :viridis,
        camera = (45, 45),
        colorbar = false,
        zlims = (vy_min, vy_max) # Control vertical axis limits
        )
    # Combine all plots into a single figure and save
    p = plot(
        p1, p2, p4, p5,
        layout = (2, 3),
        size = (1200, 800), # Increased size to give more room
        margin = 5Plots.mm, # Add margin between plots
        plot_title = "Residual Surfaces $msg"
        )
    # Check if output directory exists, if not create it
    check_dir(OUTPUT_DIR_NAME)
    savefig(p, joinpath(OUTPUT_DIR_NAME, "residual_surfaces_$(msg)_iter$(ivcycle).png"))
    cleanup()
    return nothing
end

""" Plot the convergence history of the multigrid solver.
"""
function plot_residuals(
    resx_principle::Vector{Float64},
    resy_principle::Vector{Float64},
    resz_principle::Vector{Float64},
    resc_principle::Vector{Float64},
    nvcycles::Int64,
    ivcycle::Int64;
    msg::String="MainCycle",
    ymin::Float64=-12.0,
    ymax::Float64=-2.0
)::Nothing
    iterations = 1:nvcycles

    p = plot(iterations, resc_principle, label="Continuity", color=:black, linewidth=2, reuse=false)
    plot!(p, iterations, resx_principle, label="x-Stokes", color=:blue, linewidth=2)
    plot!(p, iterations, resy_principle, label="y-Stokes", color=:green, linewidth=2)
    plot!(p, iterations, resz_principle, label="z-Stokes", color=:red, linewidth=2)
    xlims!(1, ivcycle)
    ylims!(ymin, ymax)
    xlabel!("V-cycles")
    ylabel!("log(Residuals)")
    title!("Multigrid Convergence History $msg")
    # Check if output directory exists, if not create it
    check_dir(OUTPUT_DIR_NAME)
    savefig(p, joinpath(OUTPUT_DIR_NAME, "residuals_$(msg)_iter$(ivcycle).png"))
    cleanup()
    return nothing
end

function plot_residuals_2d(
    resx_principle::Vector{Float64},
    resy_principle::Vector{Float64},
    resc_principle::Vector{Float64},
    nvcycles::Int64,
    ivcycle::Int64;
    msg::String="MainCycle",
    ymin::Float64=-15.0,
    ymax::Float64=-2.0
)::Nothing
    iterations = 1:nvcycles

    p = plot(iterations, resc_principle, label="Continuity", color=:black, linewidth=2, reuse=false)
    plot!(p, iterations, resx_principle, label="x-Stokes", color=:blue, linewidth=2)
    plot!(p, iterations, resy_principle, label="y-Stokes", color=:green, linewidth=2)
    xlims!(1, ivcycle)
    ylims!(ymin, ymax)
    xlabel!("V-cycles")
    ylabel!("log(Residuals)")
    title!("Multigrid Convergence History $msg")
    # Check if output directory exists, if not create it
    check_dir(OUTPUT_DIR_NAME)
    savefig(p, joinpath(OUTPUT_DIR_NAME, "residuals_$(msg)_iter$(ivcycle).png"))
    cleanup()
    return nothing
end

function check_dir(output_dir_name::String)
    if !isdir(output_dir_name)
        mkdir(output_dir_name)
    end
end

function cleanup()
    GC.gc()
    try
        GR.clearws()
        GR.emergencycloseall()
    catch
    end
end

end # module 