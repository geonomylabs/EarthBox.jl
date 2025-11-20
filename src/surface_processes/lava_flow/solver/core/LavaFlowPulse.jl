""" Lava flow pulse model for simulating lava flow dynamics on topography.
"""
module LavaFlowPulse

""" Extrude a lava flow pulse onto topography.

# Arguments
- `topo_gridx`: The x-coordinates (meters) of the topography grid
- `topo_gridy`: The y-coordinates (meters) of the topography grid
- `lava_thickness`: The thickness of lava on the topography grid (meters)
- `eruption_thickness`: The thickness of lava at the eruption point (meters)
- `residual_lava_thickness`: The residual thickness of lava (meters)
- `x_eruption`: The x-coordinate of the eruption point (meters)
- `tolerance`: The tolerance for the lava flow model
- `nmax`: The maximum number of iterations for the lava flow model

# Returns
- `icount`: The number of iterations for the lava flow model
"""
function lava_flow_pulse(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    lava_thickness::Vector{Float64},
    eruption_thickness::Float64,
    residual_lava_thickness::Float64,
    x_eruption::Float64;
    tolerance::Float64=1e-6,
    nmax::Int=20
)::Int
    xnum = length(topo_gridx)
    dx = topo_gridx[2] - topo_gridx[1]
    eruption_node = floor(Int, x_eruption / dx) + 1
    if eruption_node < 1 || eruption_node > xnum
        eruption_node = 1
    end

    lava_thickness[eruption_node] += eruption_thickness
    lava_thickness_old = Vector{Float64}(undef, xnum) #fill(1e38, xnum)
    Threads.@threads for i in 1:xnum
        lava_thickness_old[i] = 1e38
    end

    max_diff = calculate_max_difference(lava_thickness, lava_thickness_old)
    sorted_indices = radiate_indices(xnum, eruption_node)
    icount = 0

    while max_diff > tolerance && icount < nmax
        lava_thickness_old .= lava_thickness
        for mm in 1:xnum
            i = sorted_indices[mm]
            if 1 < i < xnum
                thickness_left = lava_thickness[i-1]
                thickness_active = lava_thickness[i]
                thickness_right = lava_thickness[i+1]

                topo_left = topo_gridy[i-1]
                topo_active = topo_gridy[i]
                topo_right = topo_gridy[i+1]

                y_surface_left = topo_left - thickness_left
                y_surface_active = topo_active - thickness_active
                y_surface_right = topo_right - thickness_right

                delta_elevation_left = max(0, y_surface_left - y_surface_active)
                delta_elevation_right = max(0, y_surface_right - y_surface_active)

                delta_elevation_nonzero_minimum = calculate_minimum_nonzero_delta_elevation(
                    delta_elevation_left, delta_elevation_right
                )

                delta_elevation_total = abs(delta_elevation_left) + abs(delta_elevation_right)

                if delta_elevation_total > 0 && thickness_active > residual_lava_thickness
                    total_potential_outflow_thickness = calculate_total_potential_outflow_thickness(
                        lava_thickness[i],
                        residual_lava_thickness,
                        delta_elevation_nonzero_minimum
                    )

                    total_outflow_thickness = 0.0
                    if delta_elevation_left > 0
                        outflow_thickness_to_left = calculate_outflow_thickness(
                            total_potential_outflow_thickness,
                            delta_elevation_left,
                            delta_elevation_total
                        )
                        total_outflow_thickness += outflow_thickness_to_left
                        lava_thickness[i-1] += outflow_thickness_to_left
                    end

                    if delta_elevation_right > 0
                        outflow_thickness_to_right = calculate_outflow_thickness(
                            total_potential_outflow_thickness,
                            delta_elevation_right,
                            delta_elevation_total
                        )
                        lava_thickness[i+1] += outflow_thickness_to_right
                        total_outflow_thickness += outflow_thickness_to_right
                    end

                    lava_thickness[i] -= total_outflow_thickness
                end
            end
        end
        max_diff = calculate_max_difference(lava_thickness, lava_thickness_old)
        icount += 1
    end
    return icount
end

""" Radiate indices from a reference index.

# Arguments
- `size`: Size of the array
- `ref_index`: Reference index to radiate from

# Returns
- Array of indices sorted by distance from reference index
"""
function radiate_indices(size::Int, ref_index::Int)::Vector{Int}
    indices = collect(1:size)
    distances = abs.(indices .- ref_index)
    return indices[sortperm(distances)]
end

""" Calculate the maximum difference between new and old lava thicknesses.

# Arguments
- `lava_thickness`: Current lava thickness array
- `lava_thickness_old`: Previous iteration's lava thickness array

# Returns
- Maximum absolute difference between arrays
"""
function calculate_max_difference(
    lava_thickness::Vector{Float64},
    lava_thickness_old::Vector{Float64}
)::Float64
    return maximum(abs.(lava_thickness .- lava_thickness_old))
end

""" Calculate thickness of lava flow out of a node.

# Arguments
- `total_potential_outflow_thickness`: Total thickness that can flow out (m)
- `delta_elevation`: Elevation difference between nodes (m)
- `delta_elevation_total`: Total elevation difference (m)
- `limit_outflow_with_elevation`: Whether to limit outflow by elevation

# Returns
- Calculated outflow thickness
"""
@inline function calculate_outflow_thickness(
    total_potential_outflow_thickness::Float64,
    delta_elevation::Float64,
    delta_elevation_total::Float64;
    limit_outflow_with_elevation::Bool=true
)::Float64
    outflow_thickness = total_potential_outflow_thickness * 
                       delta_elevation / delta_elevation_total
    
    if limit_outflow_with_elevation && outflow_thickness > delta_elevation
        outflow_thickness = delta_elevation
    end
    
    return outflow_thickness
end

""" Calculate total potential outflow thickness.

# Arguments
- `lava_thickness`: Current lava thickness
- `residual_lava_thickness`: Minimum residual thickness
- `delta_elevation_nonzero_minimum`: Minimum elevation difference
- `limit_outflow_with_elevation`: Whether to limit by elevation
- `outflow_reduction_factor`: Factor to reduce outflow

# Returns
- Total potential outflow thickness
"""
@inline function calculate_total_potential_outflow_thickness(
    lava_thickness::Float64,
    residual_lava_thickness::Float64,
    delta_elevation_nonzero_minimum::Float64;
    limit_outflow_with_elevation::Bool=true,
    outflow_reduction_factor::Float64=2/3
)::Float64
    max_potential_outflow_thickness = max(0, lava_thickness - residual_lava_thickness)

    if limit_outflow_with_elevation
        max_potential_outflow_thickness = min(
            max_potential_outflow_thickness,
            delta_elevation_nonzero_minimum
        )
    end

    return max_potential_outflow_thickness * outflow_reduction_factor
end

""" Calculate minimum nonzero delta elevation.

# Arguments
- `delta_elevation_left`: Left elevation difference
- `delta_elevation_right`: Right elevation difference

# Returns
- Minimum nonzero delta elevation
"""
@inline function calculate_minimum_nonzero_delta_elevation(
    delta_elevation_left::Float64,
    delta_elevation_right::Float64
)::Float64
    if delta_elevation_left > 0 && delta_elevation_right > 0
        return (delta_elevation_left + delta_elevation_right) / 2
    elseif delta_elevation_right > 0
        return delta_elevation_right
    else
        return delta_elevation_left
    end
end

end # module 