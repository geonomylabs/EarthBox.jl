module SedimentTransportTest

import Plots
import EarthBox.ConversionFuncs: mm_per_yr_to_meters_per_seconds as mm_yr_to_m_s
import EarthBox.ConversionFuncs: years_to_seconds
import EarthBox.ConversionFuncs: meters_per_year_to_meters_per_seconds as m_yr_to_m_s
import EarthBox.ConversionFuncs: meters_squared_per_year_to_meters_squared_per_second as m2_yr_to_m2_s
import EarthBox.ConversionFuncs: seconds_to_years
import EarthBox.DataStructures: SedimentTransportParameters
import EarthBox.SurfaceProcesses.SedimentTransport: SedimentTransportSolverManager


""" Structure to hold topography grid geometry parameters.

# Fields
- `xsize::Float64`: x-direction size of topo grid (meters)
- `dx::Float64`: spacing of topo grid nodes (meters)
- `xo_topo::Float64`: initial x-coordinate of positive topography (meters)
- `width_topo::Float64`: width of positive topography (meters)
- `height_topo::Float64`: height of positive topography (meters)
- `width_basin::Float64`: width of adjacent basin (meters)
- `depth_basin::Float64`: depth of adjacent basin (meters)
- `y_sealevel::Float64`: sea level (meters)
- `sediment_thickness_original::Float64`: original sediment thickness (meters)
"""
Base.@kwdef struct TopoGeometry
    xsize::Float64 = 500_000.0
    dx::Float64 = 500.0
    xo_topo::Float64 = 200_000.0
    width_topo::Float64 = 100_000.0
    height_topo::Float64 = 1000.0
    width_basin::Float64 = 100_000.0
    depth_basin::Float64 = 2000.0
    y_sealevel::Float64 = 0.0
    sediment_thickness_original::Float64 = 2000.0
end

""" Structure to hold topography grid data.

# Fields
- `topo_gridx::Vector{Float64}`: x-coordinate of topo grid nodes (meters)
- `topo_gridy::Vector{Float64}`: topography at grid nodes (meters)
- `gridx_b::Vector{Float64}`: x-coordinate of basic grid nodes (meters)
"""
struct TopoGrids
    topo_gridx::Vector{Float64}
    topo_gridy::Vector{Float64}
    gridx_b::Vector{Float64}
end

""" Test downhill diffusion for multiple time steps using a regular grid.
"""
function run_test()
    pelagic_sedimentation_rate = mm_yr_to_m_s(0.3) # m/s

    sediment_transport_parameters = SedimentTransportParameters(
        m2_yr_to_m2_s(0.25),
        m_yr_to_m_s(1.0), 
        1e-2,
        m2_yr_to_m2_s(1e2),
        2000.0,
        years_to_seconds(5_000.0),
        years_to_seconds(5.0*1e6),
        0.4,
        1/2000.0
    )

    topo_geometry = TopoGeometry(
        xsize=500_000.0, dx=500.0, xo_topo=200_000.0, width_topo=100_000.0,
        height_topo=1000.0, width_basin=100_000.0, depth_basin=2000.0,
        y_sealevel=0.0, sediment_thickness_original=2000.0
    )

    (
        topo_gridx, 
        topo_gridy_initial, 
        gridx_b, 
        sediment_thickness_initial
    ) = initialize_grids!(topo_geometry)

    transport_solver = SedimentTransportSolverManager.SedimentTransportSolver(
        (gridx_b[1], gridx_b[end]),
        topo_gridx,
        topo_gridy_initial,
        sediment_thickness_initial,
        sediment_transport_parameters,
        topo_geometry.y_sealevel,
        pelagic_sedimentation_rate,
        use_collections=true,
        use_print_debug=true,
        use_constant_diffusivity=false,
        use_compaction_correction=true
    )

    SedimentTransportSolverManager.run_sediment_transport_time_steps!(transport_solver)

    plot_timesteps = true
    if plot_timesteps
        make_plots(transport_solver, topo_gridx, nskip=5)
    end

    plot_displacement = true
    if plot_displacement && transport_solver.use_compaction_correction
        println(">> Min displacement: ", 
                minimum(transport_solver.compaction_displacement_max))
        println(">> Max displacement: ", 
                maximum(transport_solver.compaction_displacement_max))
        make_max_displacement_plot(
            transport_solver, topo_gridx, sediment_thickness_initial)
    end
end

function initialize_grids!(
    topo_geometry::TopoGeometry
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    xsize = topo_geometry.xsize
    dx = topo_geometry.dx
    xnum = floor(Int, xsize/dx) + 1
    topo_gridx = collect(range(0, xsize, length=xnum))
    topo_gridy = zeros(xnum)
    add_topography_to_topo_grid!(topo_gridx, topo_gridy, topo_geometry)
    add_basin_to_topo_grid!(topo_gridx, topo_gridy, topo_geometry)
    gridx_b = collect(range(0, xsize, length=xnum))
    sediment_thickness_original = 
        fill(topo_geometry.sediment_thickness_original, xnum)
    zero_out_thickness_outside_of_basin!(
        topo_gridx, sediment_thickness_original, topo_geometry)
    return topo_gridx, topo_gridy, gridx_b, sediment_thickness_original
end

function add_topography_to_topo_grid!(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    topo_geometry::TopoGeometry
)::Nothing
    xo_topo = topo_geometry.xo_topo
    width_topo = topo_geometry.width_topo
    height_topo = topo_geometry.height_topo
    
    for i in eachindex(topo_gridx)
        x_grid = topo_gridx[i]
        if xo_topo <= x_grid <= xo_topo + width_topo
            topo_gridy[i] -= height_topo
        end
    end
    return nothing
end

function add_basin_to_topo_grid!(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    topo_geometry::TopoGeometry
)::Nothing
    xo_topo = topo_geometry.xo_topo
    width_topo = topo_geometry.width_topo
    width_basin = topo_geometry.width_basin
    depth_basin = topo_geometry.depth_basin
    xo_basin = xo_topo + width_topo
    xf_basin = xo_basin + width_basin
    
    for i in eachindex(topo_gridx)
        x_grid = topo_gridx[i]
        if xo_basin <= x_grid <= xf_basin
            topo_gridy[i] += depth_basin
        end
    end
    return nothing
end

function zero_out_thickness_outside_of_basin!(
    topo_gridx::Vector{Float64},
    sediment_thickness_original::Vector{Float64},
    topo_geometry::TopoGeometry
)::Nothing
    xo_topo = topo_geometry.xo_topo
    width_topo = topo_geometry.width_topo
    width_basin = topo_geometry.width_basin
    xo_basin = xo_topo + width_topo
    xf_basin = xo_basin + width_basin
    
    for i in eachindex(topo_gridx)
        x_grid = topo_gridx[i]
        if x_grid <= xo_basin || x_grid >= xf_basin
            sediment_thickness_original[i] = 0.0
        end
    end
    return nothing
end

function make_plots(
    transport_solver::SedimentTransportSolverManager.SedimentTransportSolver,
    topo_gridx::Vector{Float64};
    nskip::Int=0
)::Nothing
    transport_timestep = transport_solver.transport_timestep
    basement_collection = transport_solver.collections.basement_collection
    topo_collection = transport_solver.collections.topo_collection
    water_depth_collection = transport_solver.collections.water_depth_collection
    divides_collection = transport_solver.collections.divides_collection

    bsmt_gridy = basement_collection[0]

    dpi = 150
    figsize = (5, 5)
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)

    for (i, topo_gridy) in topo_collection
        if nskip != 0 && i % nskip != 0
            continue
        end
        println(">> Making plot for time step: ", i)
        
        p = Plots.plot(
            topo_gridx./1000.0, topo_gridy, label="Topography", 
            color=:blue, linewidth=2.0, linestyle=:solid,
            aspect_ratio=:auto, dpi=dpi, size=figsize_pixels
            )
        Plots.plot!(
            p, topo_gridx./1000.0, bsmt_gridy, label="Initial Basement", 
            color=:red, linewidth=2.0, linestyle=:solid
            )
        
        plot_water_depth = true
        if plot_water_depth
            Plots.plot!(
                p, topo_gridx./1000.0, water_depth_collection[i], label="Water Depth", 
                color=:cyan, linewidth=2.0, linestyle=:dash
                )
        end
        
        Plots.plot!(
            p, divides_collection[i]./1000.0, zeros(length(divides_collection[i])),
            seriestype=:scatter, label="Drainage Divides", color=:green,
            marker=:circle, markersize=4.0, markerstrokecolor=:transparent, 
            markerstrokewidth=0.0
            )
        
        Plots.xlabel!("x (km)")
        Plots.ylabel!("Topography (meters)")
        time_stamp = make_time_stamp(i, transport_timestep)
        Plots.title!("Downhill Diffusion: " * time_stamp)
        Plots.ylims!(-2500, 2500)
        Plots.yaxis!(:flip)
        
        plot_name = "downhill_diffusion_" * string(i) * ".png"
        Plots.savefig(p, plot_name)
    end
    return nothing
end

function make_max_displacement_plot(
    transport_solver::SedimentTransportSolverManager.SedimentTransportSolver,
    topo_gridx::Vector{Float64},
    sediment_thickness_initial::Vector{Float64}
)::Nothing
    compaction_displacement_max = transport_solver.compaction_displacement_max
    sediment_thickness_initial_compacted = transport_solver.sediment_thickness_initial_compacted

    dpi = 150
    figsize = (5, 5)
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)

    p = Plots.plot(
        topo_gridx, sediment_thickness_initial,
        label="Initial Sediment Thickness (m)", color=:blue,
        dpi=dpi, size=figsize_pixels, linewidth=2.0, linestyle=:solid
        )
    Plots.plot!(
        p, topo_gridx, sediment_thickness_initial_compacted,
          label="Initial Sediment Thickness Compacted (m)", color=:green,
          linewidth=2.0, linestyle=:solid
          )
    Plots.plot!(
        p, topo_gridx, compaction_displacement_max,
          label="Compaction Displacement (m)", color=:red,
          linewidth=2.0, linestyle=:solid
          )
    
    Plots.xlabel!("x (m)")
    Plots.ylabel!("Displacement (m)")
    Plots.title!("Maximum Compaction Displacement in Original Layer")
    Plots.ylims!(0, 2500.0)
    
    Plots.savefig(p, "compaction_displacement.png")
    return nothing
end

function make_time_stamp(i::Int, timestep::Float64)::String
    model_time_yr = seconds_to_years(i * timestep)
    model_time_my = model_time_yr/1e6
    return string(Float64(model_time_my)) * " Myr"
end

end
