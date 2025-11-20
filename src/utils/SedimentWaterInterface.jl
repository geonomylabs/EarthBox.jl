module SedimentWaterInterface

import EarthBox.MathTools: linear_interp_at_x_location
import EarthBox.MathTools: linear_interp_bisection

""" Calculate submud depth.

# Arguments
- `x_location`: X location of marker in meters
- `y_marker`: Y location of marker in meters
- `gridt`: Topography grid array with X-coordinates stored in gridt[1, :] and
    Y-coordinates stored in grid[2, :]

# Returns
- `y_submud`: Distance below mudline (meters)
"""
function calculate_submud_depth_for_marker(
    x_location::Float64,
    y_marker::Float64,
    gridt::Matrix{Float64}
)::Float64
    y_mudline = get_depth(x_location, gridt)
    y_submud = y_marker - y_mudline
    return y_submud
end

""" Calculate y-coordinate of sediment water interface (mudline).

# Arguments
- `x_location`: X location of marker in meters
- `gridt`: Topography grid array with X-coordinates stored in gridt[1, :] and
    Y-coordinates stored in grid[2, :]

# Returns
- `y_mudline`: Y-coordinate (meters) of the sediment water interface (i.e. mudline)
"""
function get_depth(
    x_location::Float64,
    gridt::Matrix{Float64}
)::Float64
    topo_gridx, topo_gridy = get_topo_coordinates(gridt)
    y_mudline = linear_interp_at_x_location(x_location, topo_gridx, topo_gridy)
    return y_mudline
end

function get_topo_coordinates(
    gridt::Matrix{Float64}
)::Tuple{Vector{Float64},Vector{Float64}}
    toponum = size(gridt, 2)
    topo_gridx = Vector{Float64}(undef, toponum)
    topo_gridy = Vector{Float64}(undef, toponum)
    for j in 1:toponum
        topo_gridx[j] = gridt[1, j]
        topo_gridy[j] = gridt[2, j]
    end
    return topo_gridx, topo_gridy
end

""" Calculate y-coordinate of sediment water interface (mudline).

# Arguments
- `x_location`: X location of marker in meters
- `topo_gridx`: X-coordinates of topography grid (meters)
- `topo_gridy`: Y-coordinates of topography grid (meters)

# Returns
- `y_mudline`: Y-coordinate (meters) of the sediment water interface (i.e. mudline)
"""
function get_depth_opt(
    x_location::Float64,
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64}
)::Float64
    y_mudline = linear_interp_bisection(topo_gridx, topo_gridy, x_location)
    return y_mudline
end

end # module 