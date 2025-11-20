module LithostaticPressure

""" Calculate lithostatic pressure from marker swarm using a cellular mesh.

The cellular mesh is a vertical stack of cells. Density is defined in each cell by 
averaging the density of markers within the cell.

Pressure is calculated using the average density of the cell above the node and the 
gravitational acceleration.

# Arguments
- `marker_x`: Marker x-coordinates in meters
- `marker_y`: Marker y-coordinates in meters
- `marker_rho`: Marker densities in kg/m^3
- `ysize`: Size of the computation model in the y-direction in meters
- `x_location`: x-coordinate location of column lateral mid-point in meters
- `y_start`: Starting location of the column in y-direction where y increases with depth in meters
- `cell_thickness_y`: Thickness of cells used to calculate pressure by averaging marker density
- `cell_thickness_x`: Width of cells used to calculate pressure by averaging marker density

# Returns
- `gridy`: y-coordinates of the grid cell centers in meters
- `density_gridy_from_markers`: Density grid from markers in kg/m^3
- `pressure_gridy_from_markers`: Pressure grid from markers in Pascals
"""
function calculate_lithostatic_pressure_from_marker_swarm(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_rho::Vector{Float64},
    ysize::Float64,
    x_location::Float64,
    y_start::Float64,
    cell_thickness_y::Float64,
    cell_thickness_x::Float64
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    
    column_height = ysize - y_start
    ncells = floor(Int, column_height / cell_thickness_y)
    nnodes = ncells + 1

    (
        marker_x_filtered, marker_y_filtered, marker_rho_filtered
    ) = filter_markers_for_column(
        marker_x, marker_y, marker_rho, ysize,
        x_location, y_start, cell_thickness_x
    )

    gridy, density_gridy_from_markers = calculate_density_grid(
        marker_x_filtered, marker_y_filtered,
        marker_rho_filtered, y_start, cell_thickness_y, ncells
    )

    pressure_gridy_from_markers = zeros(Float64, nnodes)

    gravity_m_s_s = 9.81
    for i in 1:nnodes
        if i > 1
            density_cell_above_node = (
                density_gridy_from_markers[i-1] + 
                density_gridy_from_markers[i]
            ) / 2.0

            delta_pressure = density_cell_above_node * gravity_m_s_s * cell_thickness_y

            pressure_gridy_from_markers[i] = pressure_gridy_from_markers[i-1] + delta_pressure
        end
    end

    return gridy, density_gridy_from_markers, pressure_gridy_from_markers
end

function filter_markers_for_column(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_rho::Vector{Float64},
    ysize::Float64,
    x_location::Float64,
    y_start::Float64,
    dx::Float64
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}}
    
    nmarkers = length(marker_x)
    marker_x_tmp = zeros(Float64, nmarkers)
    marker_y_tmp = zeros(Float64, nmarkers)
    marker_rho_tmp = zeros(Float64, nmarkers)

    nmarkers_filtered = 0
    for i in 1:nmarkers
        if marker_x[i] >= x_location - dx/2 && marker_x[i] <= x_location + dx/2
            if marker_y[i] >= y_start && marker_y[i] <= ysize
                nmarkers_filtered += 1
                marker_x_tmp[nmarkers_filtered] = marker_x[i]
                marker_y_tmp[nmarkers_filtered] = marker_y[i]
                marker_rho_tmp[nmarkers_filtered] = marker_rho[i]
            end
        end
    end

    marker_x_filtered = zeros(Float64, nmarkers_filtered)
    marker_y_filtered = zeros(Float64, nmarkers_filtered)
    marker_rho_filtered = zeros(Float64, nmarkers_filtered)
    
    for i in 1:nmarkers_filtered
        marker_x_filtered[i] = marker_x_tmp[i]
        marker_y_filtered[i] = marker_y_tmp[i]
        marker_rho_filtered[i] = marker_rho_tmp[i]
    end

    return marker_x_filtered, marker_y_filtered, marker_rho_filtered
end

function calculate_density_grid(
    marker_x_filtered::Vector{Float64},
    marker_y_filtered::Vector{Float64},
    marker_rho_filtered::Vector{Float64},
    y_start::Float64,
    dy::Float64,
    ncells::Int
)::Tuple{Vector{Float64}, Vector{Float64}}
    
    gridy_cell_centers = zeros(Float64, ncells)
    density_gridy_cell_centers = zeros(Float64, ncells)

    nmarkers_filtered = length(marker_x_filtered)
    for i in 1:ncells
        y_cell_top = y_start + (i-1)*dy
        y_cell_bottom = y_cell_top + dy
        gridy_cell_centers[i] = (y_cell_top + y_cell_bottom)/2.0 - y_start

        density_sum = 0.0
        nmarkers_in_cell = 0
        for j in 1:nmarkers_filtered
            y_marker = marker_y_filtered[j]
            density = marker_rho_filtered[j]
            if y_marker >= y_cell_top && y_marker < y_cell_bottom
                density_sum += density
                nmarkers_in_cell += 1
            end
        end

        if nmarkers_in_cell > 0
            density_gridy_cell_centers[i] = density_sum / nmarkers_in_cell
        else
            error("No markers found in cell at cell index $i, cell top $y_cell_top, " *
                  "and cell bottom $y_cell_bottom")
        end
    end

    # Density on nodes
    nnodes = ncells + 1
    gridy = zeros(Float64, nnodes)
    density_gridy = zeros(Float64, nnodes)

    for i in 1:nnodes
        gridy[i] = (i-1)*dy
        if i == 1
            density_gridy[i] = density_gridy_cell_centers[1]
        elseif i == nnodes
            density_gridy[i] = density_gridy_cell_centers[ncells]
        else
            density_gridy[i] = (
                density_gridy_cell_centers[i-1] + 
                density_gridy_cell_centers[i]
            ) / 2.0
        end
    end

    return gridy, density_gridy
end

end # module 