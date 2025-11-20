module MeltDrainageDividesTest

using Plots
import EarthBox.MeltModel.Drainage: calculate_drainage_divides

struct TopoGeometry
    xsize::Float64  # x-direction size of topo grid (meters)
    dx::Float64     # spacing of topo grid nodes (meters)
    amplitude::Float64  # amplitude of topography (meters)
    wavelength::Float64  # wavelength of topography (meters)
    flat_threshold::Float64  # threshold for flat tops (meters)
end

struct TopoGrids
    topo_gridx::Vector{Float64}  # x-coordinate of topo grid nodes (meters)
    topo_gridy::Vector{Float64}  # topography at grid nodes (meters)
end

function run_test(;use_plotter::Bool=false)::Vector{Float64}
    topo_geometry = TopoGeometry(
        500_000.0,  # xsize
        500.0,      # dx
        5_000.0,    # amplitude
        100_000.0,  # wavelength
        2_500.0     # flat_threshold
    )
    topo_grids = initialize_grids(topo_geometry)
    divides_x = calculate_drainage_divides(topo_grids.topo_gridx, topo_grids.topo_gridy)
    if use_plotter
        make_plots(topo_grids.topo_gridx, topo_grids.topo_gridy, divides_x)
    end
    return divides_x
end

function initialize_grids(topo_geometry::TopoGeometry)::TopoGrids
    # Unpack parameters
    xsize = topo_geometry.xsize
    dx = topo_geometry.dx
    xnum = floor(Int, xsize/dx) + 1
    topo_gridx = collect(range(0, xsize, length=xnum))
    topo_gridy = zeros(xnum)
    add_topography_to_topo_grid!(topo_gridx, topo_gridy, topo_geometry)
    return TopoGrids(topo_gridx, topo_gridy)
end

function add_topography_to_topo_grid!(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    topo_geometry::TopoGeometry
)::Nothing
    amplitude = topo_geometry.amplitude
    wavelength = topo_geometry.wavelength
    flat_threshold = topo_geometry.flat_threshold

    sine_values = amplitude .* sin.(topo_gridx ./ wavelength) .+
                 amplitude .* sin.(3π .* topo_gridx ./ wavelength)

    # Flatten the bottoms
    sine_values[sine_values .> flat_threshold] .= flat_threshold
    # Flatten tops
    sine_values[sine_values .< -flat_threshold * 2.0] .= -flat_threshold * 2.0
    # Update topo_gridy
    topo_gridy .= sine_values
    return nothing
end

function make_plots(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    divides_x::Vector{Float64}
)::Nothing
    p = plot(
        topo_gridx, topo_gridy,
        label="Topography",
        color=:blue,
        xlabel="x (m)",
        ylabel="Topography",
        title="Melt drainage divides test",
        ylims=(-10000, 10000),
        yflip=true
    )
    
    scatter!(
        p,
        divides_x,
        zeros(length(divides_x)),
        label="Drainage Divides",
        color=:green,
        marker=:circle
    )
    
    savefig(p, "melt_drainage_divides_test.png")
    return nothing
end

end # module 