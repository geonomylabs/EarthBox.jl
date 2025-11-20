"""
    ResidualPlottingCM

Module for plotting residual surfaces for Stokes and continuity equations in 3D.
This module provides functions to create surface plots similar to the MATLAB 
implementation in Variable_viscosity3D_Multigrid.m using CairoMakie instead of Plots.
"""
module ResidualPlottingCM

using CairoMakie
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
    msg::String="PrincipleLevel1"
)::Nothing
    min_delta = 1e-15
    # Calculate the middle z-slice index
    znum = level_data.grid.parameters.geometry.znum.value
    vy = level_data.vy.array

    z_mid = div(znum - 1, 2) + 1
    
    # Create coordinate grids for each residual array separately
    # For CairoMakie, we need to ensure coordinate arrays match data dimensions exactly
    resx0_min = minimum(resx0)
    resx0_max = maximum(resx0)
    delta = resx0_max - resx0_min
    #println("resx0_min = $resx0_min, resx0_max = $resx0_max, delta = $delta")
    if delta < min_delta
        resx0_max = resx0_min + min_delta
    end
    
    resy0_min = minimum(resy0)
    resy0_max = maximum(resy0)
    delta = resy0_max - resy0_min
    #println("resy0_min = $resy0_min, resy0_max = $resy0_max, delta = $delta")
    if delta < min_delta
        resy0_max = resy0_min + min_delta
    end
    
    resz0_min = minimum(resz0)
    resz0_max = maximum(resz0)
    delta = resz0_max - resz0_min
    #println("resz0_min = $resz0_min, resz0_max = $resz0_max, delta = $delta")
    if delta < min_delta
        resz0_max = resz0_min + min_delta
    end
    
    resc0_min = minimum(resc0)
    resc0_max = maximum(resc0)
    delta = resc0_max - resc0_min
    #println("resc0_min = $resc0_min, resc0_max = $resc0_max, delta = $delta")
    if delta < min_delta
        resc0_max = resc0_min + min_delta
    end

    vy_min = minimum(vy)
    vy_max = maximum(vy)
    delta = vy_max - vy_min
    #println("vy_min = $vy_min, vy_max = $vy_max, delta = $delta")
    if delta < min_delta
        vy_max = vy_min + min_delta
    end
    
    # Create figure with CairoMakie
    fig = Figure(size = (1200, 800))
    
    # Create subplots
    ax1 = Axis3(
        fig[1, 1], 
        title = ivcycle > 0 ? "x-Stokes residual, cycle = $ivcycle" : "x-Stokes residual",
        xlabel = "x", ylabel = "y", zlabel = "residual x-Stokes",
        limits = ((1, size(resx0, 2)), (1, size(resx0, 1)), (resx0_min, resx0_max))
        )
    surface!(ax1, resx0[:, :, z_mid], colormap = :viridis, colorrange = (resx0_min, resx0_max))
    
    ax2 = Axis3(
        fig[1, 2], 
        title = ivcycle > 0 ? "y-Stokes residual, cycle = $ivcycle" : "y-Stokes residual",
        xlabel = "x", ylabel = "y", zlabel = "residual y-Stokes",
        limits = ((1, size(resy0, 2)), (1, size(resy0, 1)), (resy0_min, resy0_max))
        )
    surface!(ax2, resy0[:, :, z_mid], colormap = :viridis, colorrange = (resy0_min, resy0_max))
    
    ax3 = Axis3(
        fig[1, 3], 
        title = ivcycle > 0 ? "z-Stokes residual, cycle = $ivcycle" : "z-Stokes residual",
        xlabel = "x", ylabel = "y", zlabel = "residual z-Stokes",
        limits = ((1, size(resz0, 2)), (1, size(resz0, 1)), (resz0_min, resz0_max))
        )
    surface!(ax3, resz0[:, :, z_mid], colormap = :viridis, colorrange = (resz0_min, resz0_max))
    
    ax4 = Axis3(
        fig[2, 1], 
        title = ivcycle > 0 ? "Continuity residual, cycle = $ivcycle" : "Continuity residual",
        xlabel = "x", ylabel = "y", zlabel = "residual continuity",
        limits = ((1, size(resc0, 2)), (1, size(resc0, 1)), (resc0_min, resc0_max))
        )
    surface!(ax4, resc0[:, :, z_mid], colormap = :viridis, colorrange = (resc0_min, resc0_max))
    
    ax5 = Axis3(
        fig[2, 2], 
        title = "Vy at i = $z_mid (z_mid)",
        xlabel = "x", ylabel = "y", zlabel = "Vy (m/s)",
        limits = ((1, size(vy, 2)), (1, size(vy, 1)), (vy_min, vy_max))
        )
    surface!(ax5, vy[:, :, z_mid], colormap = :viridis, colorrange = (vy_min, vy_max))
    
    # Add overall title
    Label(fig[0, :], "Residual Surfaces $msg", fontsize = 16)
    
    # Check if output directory exists, if not create it
    check_dir(OUTPUT_DIR_NAME)
    save(joinpath(OUTPUT_DIR_NAME, "residual_surfaces_iter$(ivcycle).png"), fig)
    return nothing
end

function plot_residual_surfaces(
    resx0::Array{Float64, 2},
    resy0::Array{Float64, 2}, 
    resc0::Array{Float64, 2},
    level_data::LevelData2d,
    ivcycle::Int64=0;
    msg::String="PrincipleLevel1"
)::Nothing
    min_delta = 1e-15

    # Calculate the middle z-slice index
    vy = level_data.vy.array

    # Create coordinate grids for each residual array separately
    # For CairoMakie, we need to ensure coordinate arrays match data dimensions exactly
    resx0_min = minimum(resx0)
    resx0_max = maximum(resx0)
    delta = resx0_max - resx0_min
    #println("resx0_min = $resx0_min, resx0_max = $resx0_max, delta = $delta")
    if delta < min_delta
        resx0_max = resx0_min + min_delta
    end
    
    resy0_min = minimum(resy0)
    resy0_max = maximum(resy0)
    delta = resy0_max - resy0_min
    #println("resy0_min = $resy0_min, resy0_max = $resy0_max, delta = $delta")
    if delta < min_delta
        resy0_max = resy0_min + min_delta
    end
    
    resc0_min = minimum(resc0)
    resc0_max = maximum(resc0)
    delta = resc0_max - resc0_min
    #println("resc0_min = $resc0_min, resc0_max = $resc0_max, delta = $delta")
    if delta < min_delta
        resc0_max = resc0_min + min_delta
    end

    vy_min = minimum(vy)
    vy_max = maximum(vy)
    delta = vy_max - vy_min
    #println("vy_min = $vy_min, vy_max = $vy_max, delta = $delta")
    if delta < min_delta
        vy_max = vy_min + min_delta
    end
    
    # Create figure with CairoMakie
    fig = Figure(size = (1200, 800))
    
    # Create subplots
    ax1 = Axis3(
        fig[1, 1], 
        title = ivcycle > 0 ? "x-Stokes residual, cycle = $ivcycle" : "x-Stokes residual",
        xlabel = "x", ylabel = "y", zlabel = "residual x-Stokes",
        limits = ((1, size(resx0, 2)), (1, size(resx0, 1)), (resx0_min, resx0_max))
        )
    surface!(ax1, resx0, colormap = :viridis, colorrange = (resx0_min, resx0_max))
    
    ax2 = Axis3(
        fig[1, 2], 
        title = ivcycle > 0 ? "y-Stokes residual, cycle = $ivcycle" : "y-Stokes residual",
        xlabel = "x", ylabel = "y", zlabel = "residual y-Stokes",
        limits = ((1, size(resy0, 2)), (1, size(resy0, 1)), (resy0_min, resy0_max))
        )
    surface!(ax2, resy0, colormap = :viridis, colorrange = (resy0_min, resy0_max))
    
    ax4 = Axis3(
        fig[2, 1], 
        title = ivcycle > 0 ? "Continuity residual, cycle = $ivcycle" : "Continuity residual",
        xlabel = "x", ylabel = "y", zlabel = "residual continuity",
        limits = ((1, size(resc0, 2)), (1, size(resc0, 1)), (resc0_min, resc0_max))
        )
    surface!(ax4, resc0, colormap = :viridis, colorrange = (resc0_min, resc0_max))
    
    ax5 = Axis3(
        fig[2, 2], 
        title = "Vy cycle = $ivcycle",
        xlabel = "x", ylabel = "y", zlabel = "Vy (m/s)",
        limits = ((1, size(vy, 2)), (1, size(vy, 1)), (vy_min, vy_max))
        )
    surface!(ax5, vy, colormap = :viridis, colorrange = (vy_min, vy_max))
    
    # Add overall title
    Label(fig[0, :], "Residual Surfaces $msg", fontsize = 16)
    
    # Check if output directory exists, if not create it
    check_dir(OUTPUT_DIR_NAME)
    save(joinpath(OUTPUT_DIR_NAME, "residual_surfaces_$(msg)_iter$(ivcycle).png"), fig)
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
    msg::String="PrincipleLevel1",
    ymin::Float64=-12.0,
    ymax::Float64=-2.0,
    xlabel::String="V-cycles"
)::Nothing
    iterations = 1:nvcycles

    # Check for valid data
    if any(isnan, resc_principle) || any(isnan, resx_principle) || any(isnan, resy_principle) || any(isnan, resz_principle)
        println("Warning: NaN values detected in residual data, skipping plot")
        return nothing
    end

    fig = Figure(size = (800, 600))
    ax = Axis(
        fig[1, 1], 
        xlabel = xlabel, 
        ylabel = "log(Residuals)",
        title = "Multigrid Convergence History $msg"
        )
    
    # Set limits after creating the axis to avoid text bounding box issues
    # Ensure we have valid limits (different min and max values)
    x_min = 1
    x_max = max(ivcycle, 2)  # Ensure at least 2 for valid range
    y_min = ymin
    y_max = ymax
    
    xlims!(ax, x_min, x_max)
    ylims!(ax, y_min, y_max)
    
    lines!(ax, iterations, resc_principle, color = :black, linewidth = 2, label = "Continuity")
    lines!(ax, iterations, resx_principle, color = :blue, linewidth = 2, label = "x-Stokes")
    lines!(ax, iterations, resy_principle, color = :green, linewidth = 2, label = "y-Stokes")
    lines!(ax, iterations, resz_principle, color = :red, linewidth = 2, label = "z-Stokes")
    
    axislegend(ax)
    
    # Check if output directory exists, if not create it
    check_dir(OUTPUT_DIR_NAME)
    save(joinpath(OUTPUT_DIR_NAME, "residuals_$(msg)_iter$(ivcycle).png"), fig)
    return nothing
end

function plot_residuals_2d(
    resx_principle::Vector{Float64},
    resy_principle::Vector{Float64},
    resc_principle::Vector{Float64},
    nvcycles::Int64,
    ivcycle::Int64;
    msg::String="PrincipleLevel1",
    ymin::Float64=-15.0,
    ymax::Float64=-2.0,
    xlabel::String="V-cycles"
)::Nothing
    iterations = 1:nvcycles

    # Check for valid data
    if any(isnan, resc_principle) || any(isnan, resx_principle) || any(isnan, resy_principle)
        println("Warning: NaN values detected in residual data, skipping plot")
        return nothing
    end

    fig = Figure(size = (800, 600))
    ax = Axis(
        fig[1, 1], 
        xlabel = xlabel, 
        ylabel = "log(Residuals)",
        title = "Multigrid Convergence History $msg"
        )
    
    # Set limits after creating the axis to avoid text bounding box issues
    # Ensure we have valid limits (different min and max values)
    x_min = 1
    x_max = max(ivcycle, 2)  # Ensure at least 2 for valid range
    y_min = ymin
    y_max = ymax
    
    xlims!(ax, x_min, x_max)
    ylims!(ax, y_min, y_max)
    
    lines!(ax, iterations, resc_principle, color = :black, linewidth = 2, label = "Continuity")
    lines!(ax, iterations, resx_principle, color = :blue, linewidth = 2, label = "x-Stokes")
    lines!(ax, iterations, resy_principle, color = :green, linewidth = 2, label = "y-Stokes")
    
    axislegend(ax)
    
    # Check if output directory exists, if not create it
    check_dir(OUTPUT_DIR_NAME)
    save(joinpath(OUTPUT_DIR_NAME, "residuals_$(msg)_iter$(ivcycle).png"), fig)
    return nothing
end

function check_dir(output_dir_name::String)
    if !isdir(output_dir_name)
        mkpath(output_dir_name)
    end
end

end # module 