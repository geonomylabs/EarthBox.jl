module TestGridGravity

import EarthBox.Gravity: gravity_anomaly_loop_left_edge, calculate_free_air_gravity
import Plots

mutable struct GridGeometry
    xsize::Float64
    ysize::Float64
    dx::Float64
    dy::Float64

    function GridGeometry(xsize::Float64, ysize::Float64, dx::Float64, dy::Float64)
        new(xsize, ysize, dx, dy)
    end
end

mutable struct Topography
    xo_topo::Float64
    width_topo::Float64
    height_topo::Float64
    root_thickness::Float64

    function Topography(
        xo_topo::Float64, width_topo::Float64, 
        height_topo::Float64, root_thickness::Float64
    )
        new(xo_topo, width_topo, height_topo, root_thickness)
    end
end

mutable struct Thickness
    thick_crust::Float64
    thick_air::Float64

    function Thickness(thick_crust::Float64, thick_air::Float64)
        new(thick_crust, thick_air)
    end
end

mutable struct Density
    rho_crust::Float64
    rho_mantle::Float64
    rho_air::Float64
    rho_water::Float64

    function Density(
        rho_crust::Float64, rho_mantle::Float64, 
        rho_air::Float64, rho_water::Float64
    )
        new(rho_crust, rho_mantle, rho_air, rho_water)
    end
end

function run_test()::Nothing
    (grid_geometry, thickness, density, topography) = get_inputs()
    (bgridx, bgridy, xsteps, ysteps, rho_grid, topo_gridx, topo_gridy) = 
        build_grids(grid_geometry, thickness, density, topography)
    gravity_constant = 6.6732e-11
    y_sealevel = thickness.thick_air

    t1 = time()
    gravity_anomaly_mgal = gravity_anomaly_loop_left_edge(
        y_sealevel, bgridx, bgridy, xsteps, ysteps, rho_grid,
        topo_gridx, gravity_constant
    )
    gravity_anomaly_free_air_mgal = calculate_free_air_gravity(
        y_sealevel, gravity_constant, density.rho_crust,
        topo_gridx, topo_gridy, gravity_anomaly_mgal
    )
    t2 = time()
    println(">> Calculated gravity anomaly in ", round(t2-t1, digits=4), 
        " seconds.")

    make_plots(topo_gridx, gravity_anomaly_mgal, gravity_anomaly_free_air_mgal, 
        rho_grid)
end

""" Get inputs from Fowler 2001, Figure 5.6a (from Bott, 1982).

# Returns
- `grid_geometry::GridGeometry`: Grid geometry information
- `thickness::Thickness`: Thickness information
- `density::Density`: Density information
- `topography::Topography`: Topography information
"""
function get_inputs()::Tuple{GridGeometry, Thickness, Density, Topography}
    grid_geometry = GridGeometry(
        700_000.0,
        210_000.0,
        1000.0,
        1000.0
    )

    thickness = Thickness(
        30_000.0,
        10_000.0
    )

    density = Density(
        2850.0,
        3300.0,
        0.0,
        1000.0
    )

    height_topo = 3000.0 # meters

    root_thickness = calculate_mountain_root_thickness(height_topo, density)
    println(">> root_thickness: ", root_thickness)
    topography = Topography(
        100_000.0,
        500_000.0,
        height_topo,
        root_thickness
    )

    return grid_geometry, thickness, density, topography
end

""" Build grids.

# Arguments
- `grid_geometry::GridGeometry`: Grid geometry information
- `thickness::Thickness`: Thickness information
- `density::Density`: Density information
- `topography::Topography`: Topography information

# Returns
- `bgridx::Vector{Float64}`: x-coordinates of the grid (meters)
- `bgridy::Vector{Float64}`: y-coordinates of the grid (meters)
- `xsteps::Vector{Float64}`: x step sizes (meters)
- `ysteps::Vector{Float64}`: y step sizes (meters)
- `rho_grid::Matrix{Float64}`: density grid (kg/m³)
- `topo_gridx::Vector{Float64}`: x-coordinates of topography grid (meters)
- `topo_gridy::Vector{Float64}`: y-coordinates of topography grid (meters)
"""
function build_grids(
    grid_geometry::GridGeometry,
    thickness::Thickness,
    density::Density,
    topography::Topography
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, 
    Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    
    (bgridx, bgridy, xsteps, ysteps, rho_grid, topo_gridx, topo_gridy) = 
        initialize_grids(grid_geometry, density, thickness)

    add_topography_to_topo_grid(topo_gridx, topo_gridy, topography)
    define_initial_air_layer(bgridx, bgridy, rho_grid, thickness, density)
    define_initial_crust_layer(bgridx, bgridy, rho_grid, thickness, density)
    add_topography_and_root_to_density_grid(
        bgridx, bgridy, rho_grid, thickness, density, topography)
    return bgridx, bgridy, xsteps, ysteps, rho_grid, topo_gridx, topo_gridy
end

""" Make plots.

# Arguments
- `topo_gridx::Vector{Float64}`: x-coordinates of topography grid (meters)
- `gravity_anomaly_mgal::Vector{Float64}`: Bouguer gravity anomaly (mgal)
- `gravity_anomaly_free_air_mgal::Vector{Float64}`: Free air gravity anomaly (mgal)
- `rho_grid::Matrix{Float64}`: density grid (kg/m³)
"""
function make_plots(
    topo_gridx::Vector{Float64},
    gravity_anomaly_mgal::Vector{Float64},
    gravity_anomaly_free_air_mgal::Vector{Float64},
    rho_grid::Matrix{Float64}
)::Nothing
    println(">> Plotting....")
    plot_curve = true
    if plot_curve
        p1 = Plots.plot(topo_gridx, gravity_anomaly_mgal, label="Bouguer")
        Plots.plot!(p1, topo_gridx, gravity_anomaly_free_air_mgal, 
            label="Free Air")
        Plots.xlabel!(p1, "x (m)")
        Plots.ylabel!(p1, "Gravity Anomaly (mgal)")
        Plots.title!(p1, "Test gravity grid loop.")
        #Plots.legend!(p1)
        Plots.savefig(p1, "gravity_anomaly_mgal.png")
        Plots.display(p1)
    end

    plot_grid = true
    if plot_grid
        p2 = Plots.heatmap(rho_grid, color=:viridis, aspect_ratio=:equal, yflip=true)
        Plots.xlabel!(p2, "X")
        Plots.ylabel!(p2, "Y")
        Plots.title!(p2, "2D Density Grid")
        Plots.savefig(p2, "density_grid.png")
        Plots.display(p2)
    end

    println(">> Done.")
end

""" Calculate the root thickness of the mountain.

# Arguments
- `height_topo::Float64`: Height of the mountain (meters)
- `density::Density`: Density information (kg/m³)

# Returns
- `root_thickness::Float64`: Root thickness of the mountain (meters)
"""
function calculate_mountain_root_thickness(
    height_topo::Float64,
    density::Density
)::Float64
    root_thickness = (
        height_topo * density.rho_crust / (density.rho_mantle - density.rho_crust)
    )
    return root_thickness
end

""" Make grids including coordinates, density and topography.

Density grid is initialized using mantle density.
Topography grid is initialized using the y-depth of sealevel.

# Arguments
- `grid_geometry::GridGeometry`: Grid geometry information
- `density::Density`: Density information (kg/m³)
- `thickness::Thickness`: Thickness information (meters)

# Returns
- `bgridx::Vector{Float64}`: x-coordinates of the grid (meters)
- `bgridy::Vector{Float64}`: y-coordinates of the grid (meters)
- `xsteps::Vector{Float64}`: x step sizes (meters)
- `ysteps::Vector{Float64}`: y step sizes (meters)
- `rho_grid::Matrix{Float64}`: density grid (kg/m³)
- `topo_gridx::Vector{Float64}`: x-coordinates of topography grid (meters)
- `topo_gridy::Vector{Float64}`: y-coordinates of topography grid (meters)
"""
function initialize_grids(
    grid_geometry::GridGeometry,
    density::Density,
    thickness::Thickness
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}, 
    Matrix{Float64}, Vector{Float64}, Vector{Float64}}
    
    # Unpack parameters
    xsize = grid_geometry.xsize
    ysize = grid_geometry.ysize
    dx = grid_geometry.dx
    dy = grid_geometry.dy
    rho_mantle = density.rho_mantle
    y_sealevel = thickness.thick_air

    xnum = floor(Int, xsize/dx) + 1
    ynum = floor(Int, ysize/dy) + 1

    bgridx = collect(range(0.0, xsize, length=xnum))
    bgridy = collect(range(0.0, ysize, length=ynum))

    xsteps = diff(bgridx)
    ysteps = diff(bgridy)
    # print size of xsteps and ysteps
    rho_grid = ones(Float64, ynum, xnum) * rho_mantle

    topo_gridx = collect(range(0.0, xsize, length=xnum))
    topo_gridy = ones(Float64, xnum) * y_sealevel

    return bgridx, bgridy, xsteps, ysteps, rho_grid, topo_gridx, topo_gridy
end

""" Define air layer in the density grid.

# Arguments
- `bgridx::Vector{Float64}`: x-coordinates of grid nodes (meters)
- `bgridy::Vector{Float64}`: y-coordinates of grid nodes (meters)
- `rho_grid::Matrix{Float64}`: Density at grid nodes (kg/m³) - modified in place
- `thickness::Thickness`: Thickness information
- `density::Density`: Density information
"""
function define_initial_air_layer(
    bgridx::Vector{Float64},
    bgridy::Vector{Float64},
    rho_grid::Matrix{Float64},
    thickness::Thickness,
    density::Density
)::Nothing
    ynum = length(bgridy)
    xnum = length(bgridx)
    y_sealevel = thickness.thick_air
    rho_air = density.rho_air
    
    # Julia is column-major, so j (x) is outer loop, i (y) is inner loop
    for j in 1:xnum
        for i in 1:ynum
            if bgridy[i] < y_sealevel
                rho_grid[i, j] = rho_air
            end
        end
    end
end

""" Define initial crust layer in the density grid.

# Arguments
- `bgridx::Vector{Float64}`: x-coordinates of grid nodes (meters)
- `bgridy::Vector{Float64}`: y-coordinates of grid nodes (meters)
- `rho_grid::Matrix{Float64}`: Density at grid nodes (kg/m³) - modified in place
- `thickness::Thickness`: Thickness information
- `density::Density`: Density information
"""
function define_initial_crust_layer(
    bgridx::Vector{Float64},
    bgridy::Vector{Float64},
    rho_grid::Matrix{Float64},
    thickness::Thickness,
    density::Density
)::Nothing
    ynum = length(bgridy)
    xnum = length(bgridx)
    rho_crust = density.rho_crust
    y_sealevel = thickness.thick_air
    thickness_crust = thickness.thick_crust
    y_moho = y_sealevel + thickness_crust
    
    # Julia is column-major, so j (x) is outer loop, i (y) is inner loop
    for j in 1:xnum
        for i in 1:ynum
            y_grid = bgridy[i]
            if y_sealevel <= y_grid <= y_moho
                rho_grid[i, j] = rho_crust
            end
        end
    end
end

""" Add mountain topography and root in the density grid.

# Arguments
- `bgridx::Vector{Float64}`: x-coordinates of grid nodes (meters)
- `bgridy::Vector{Float64}`: y-coordinates of grid nodes (meters)
- `rho_grid::Matrix{Float64}`: Density at grid nodes (kg/m³) - modified in place
- `thickness::Thickness`: Thickness information
- `density::Density`: Density information
- `topography::Topography`: Topography information
"""
function add_topography_and_root_to_density_grid(
    bgridx::Vector{Float64},
    bgridy::Vector{Float64},
    rho_grid::Matrix{Float64},
    thickness::Thickness,
    density::Density,
    topography::Topography
)::Nothing
    ynum = length(bgridy)
    xnum = length(bgridx)
    rho_crust = density.rho_crust
    y_sealevel = thickness.thick_air
    xo_topo = topography.xo_topo
    width_topo = topography.width_topo

    height_topo = topography.height_topo
    root_thickness = topography.root_thickness
    crust_thickness = thickness.thick_crust

    y_topo = y_sealevel - height_topo
    y_root = y_sealevel + crust_thickness + root_thickness

    # Julia is column-major, so j (x) is outer loop, i (y) is inner loop
    for j in 1:xnum
        for i in 1:ynum
            y_grid = bgridy[i]
            x_grid = bgridx[j]
            if xo_topo <= x_grid <= xo_topo + width_topo
                if y_topo <= y_grid <= y_root
                    rho_grid[i, j] = rho_crust
                end
            end
        end
    end
end

""" Add mountain topography to the topography grid.

# Arguments
- `topo_gridx::Vector{Float64}`: x-coordinates of topography grid (meters)
- `topo_gridy::Vector{Float64}`: Topography at grid nodes (meters) - modified in place
- `topography::Topography`: Topography information
"""
function add_topography_to_topo_grid(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    topography::Topography
)::Nothing
    xnum = length(topo_gridx)
    xo_topo = topography.xo_topo
    width_topo = topography.width_topo
    height_topo = topography.height_topo
    
    for i in 1:xnum
        x_grid = topo_gridx[i]
        if xo_topo <= x_grid <= xo_topo + width_topo
            topo_gridy[i] = topo_gridy[i] - height_topo
        end
    end
end

end # module
