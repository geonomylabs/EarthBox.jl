module Scatter

using Printf
using LinearAlgebra
using Statistics
using CairoMakie

import ..Filters: apply_dimensional_filter

# Type aliases for better readability
const AxesType = CairoMakie.Axis

""" Plot scatter.

The custom color map is used if it is not nothing. Otherwise, the color map
is obtained using the default colormap.

# Arguments

- axes::AxesType: 
    - Axes object.
- color_map::Union{Symbol, Any}:
    - Color map.
- x_array::Vector{Float64}: 
    - X-coordinates.
- y_array::Vector{Float64}: 
    - Y-coordinates.
- color_array::Vector{Float64}: 
    - Scalar values used to define the color of the markers.
- marker_size::Float64:
    - Marker size.
- min_value::Float64:
    - Minimum scalar value.
- max_value::Float64:
    - Maximum scalar value.
- ticks::Vector{Float64}:
    - Ticks for color bar.

# Keyword Arguments

- label::Union{String, Nothing}:
    - Label for color bar.
- custom_cmap::Union{Any, Nothing}:
    - Custom color map.
- custom_labels::Union{Vector{String}, Nothing}:
    - Custom labels for color bar.
- order_number::Int: 
    - Order number for color bar for cases with multiple color bars.
- decimation_factor::Int: 
    - Decimation factor.
- colorbar_labels_fontsize::Int:
    - Font size for color bar labels.
- colorbar_ticks_fontsize::Int:
    - Font size for color bar ticks.
- plot_dimensions::Tuple{Float64, Float64, Float64, Float64}: 
    - Plot dimensions: (xmin, xmax, ymin, ymax).
- apply_domain_filter::Bool:
    - Apply domain filter.
"""
function plot_scatter(
    axes::AxesType,
    color_map::Union{Symbol, Any},
    x_array::Vector{Float64},
    y_array::Vector{Float64},
    color_array::Vector{Float64},
    marker_size::Float64,
    min_value::Float64,
    max_value::Float64,
    ticks::Vector{Float64};
    label::Union{String, Nothing} = nothing,
    custom_cmap::Union{Any, Nothing} = nothing,
    custom_labels::Union{Vector{String}, Nothing} = nothing,
    order_number::Int = 2,
    colorbar_labels_fontsize::Int = 6,
    colorbar_ticks_fontsize::Int = 6,
    decimation_factor::Int = 1,
    plot_dimensions::Tuple{Float64, Float64, Float64, Float64} = (0.0, 1e39, 0.0, 1e39),
    apply_domain_filter::Bool = false
)::Nothing

    # Use custom colormap if provided, otherwise use default
    cmap = isnothing(custom_cmap) ? color_map : custom_cmap

    # Apply decimation if requested
    if decimation_factor > 1
        @printf(">> Decimating markers using a decimation factor of %d\n", decimation_factor)
        @printf("   >> Number of markers before decimation: %d\n", length(x_array))
        indices = 1:decimation_factor:length(x_array)
        x_array = x_array[indices]
        y_array = y_array[indices]
        color_array = color_array[indices]
        @printf("   >> Number of markers after decimation: %d\n", length(x_array))
    end

    color_array = convert(Vector{Float64}, color_array)

    if apply_domain_filter
        x_array, y_array, color_array = apply_dimensional_filter(
            x_array, y_array, color_array, plot_dimensions)
    end

    scatter_plot = CairoMakie.scatter!(
        axes, 
        x_array, 
        y_array, 
        color=color_array, 
        colormap=cmap, 
        markersize=marker_size,
        marker=:circle,
        strokewidth=0.0,
        strokecolor=:black,
        colorrange=(min_value, max_value),
        rasterize=true
    )

    fig = axes.parent
    n_rows = length(fig.layout.rowsizes)
    color_bar = CairoMakie.Colorbar(fig[n_rows, order_number], scatter_plot)

    if !isnothing(label) && isnothing(custom_labels)
        color_bar.ticks = ticks
        color_bar.label = label
        color_bar.labelsize = colorbar_labels_fontsize
    end

    if !isnothing(custom_labels)
        color_bar.ticks = (ticks, custom_labels)
    end

    color_bar.ticklabelsize = colorbar_ticks_fontsize

    return nothing

end

""" Create regular grids from scatter data for contouring.

Inputs
------
- x_array_original::Vector{Float64}: Original x-coordinates
- y_array_original::Vector{Float64}: Original y-coordinates  
- values_original::Vector{Float64}: Original scalar values

Keyword Arguments
----------------
- nx::Int: Number of x grid points (default: 50)
- ny::Int: Number of y grid points (default: 50)
- decimation_factor::Int: Factor to decimate original data (default: 2)

Returns
-------
- Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}: 
  (grid_x, grid_y, grid_z) regular grids
"""
function make_grids_from_scatter(
    x_array_original::Vector{Float64},
    y_array_original::Vector{Float64},
    values_original::Vector{Float64};
    nx::Int = 50,
    ny::Int = 50,
    decimation_factor::Int = 2
)
    # Apply decimation
    indices = 1:decimation_factor:length(x_array_original)
    x_array = x_array_original[indices]
    y_array = y_array_original[indices]
    values = values_original[indices]

    # Get bounds
    x_min, x_max = extrema(x_array)
    y_min, y_max = extrema(y_array)

    # Create regular grid
    grid_x = collect(range(x_min, x_max, length=nx))
    grid_y = collect(range(y_min, y_max, length=ny))

    # Interpolate values to grid
    grid_z = interpolate_scatter_to_grid(x_array, y_array, values, grid_x, grid_y)

    # Note: Triangulation-based masking is simplified for CairoMakie
    # In production, implement proper triangulation-based masking
    #return grid_x_mat, grid_y_mat, grid_z
    return grid_x, grid_y, grid_z
end

""" Interpolate scatter data to regular grid using nearest neighbor method.

Returns
-------
- Matrix{Float64}: Interpolated values on regular grid
"""
function interpolate_scatter_to_grid(
    x_array::Vector{Float64},
    y_array::Vector{Float64}, 
    values::Vector{Float64},
    grid_x::Vector{Float64},
    grid_y::Vector{Float64}
)
    nx, ny = length(grid_x), length(grid_y)
    grid_z = zeros(Float64, ny, nx)
    
    for j in 1:ny  # Outer loop for y (second index)
        for i in 1:nx  # Inner loop for x (first index)
            # Find nearest neighbor
            distances = sqrt.((x_array .- grid_x[i]).^2 .+ (y_array .- grid_y[j]).^2)
            nearest_idx = argmin(distances)
            grid_z[j, i] = values[nearest_idx]
        end
    end
    
    return grid_z
end

"""
    point_in_polygon(x::Float64, y::Float64, poly::Matrix{Float64})

Check if point is in polygon using ray casting algorithm.

Inputs
------
- x::Float64: X-coordinate of point
- y::Float64: Y-coordinate of point  
- poly::Matrix{Float64}: Polygon vertices as 2×n matrix

Returns
-------
- Bool: true if point is inside polygon, false otherwise
"""
function point_in_polygon(x::Float64, y::Float64, poly::Matrix{Float64})
    xinters = -1e39
    n = size(poly, 2)
    inside = false
    p1x, p1y = poly[1, 1], poly[2, 1]
    
    for i in 1:n+1
        p2x, p2y = poly[1, (i-1) % n + 1], poly[2, (i-1) % n + 1]
        if y > min(p1y, p2y)
            if y <= max(p1y, p2y)
                if x <= max(p1x, p2x)
                    if p1y != p2y
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    end
                    if p1x == p2x || x <= xinters
                        inside = !inside
                    end
                end
            end
        end
        p1x, p1y = p2x, p2y
    end
    
    return inside
end

end # module
