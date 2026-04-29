""" Lava flow pulse model for simulating lava flow dynamics on topography.
"""
module LavaFlowPulse

""" Extrude a lava flow pulse onto topography.

# Arguments
- `topo_gridx`: The x-coordinates (meters) of the topography grid
- `topo_gridy`: The y-coordinates (meters) of the topography grid
- `lava_thickness`: The thickness of lava on the topography grid (meters)
- `lava_thickness_old`: Caller-provided scratch buffer of length `xnum`.
    Contents on entry are irrelevant — `fill!`-ed with `1e38` here as a
    sentinel to force the convergence loop to enter at least once.
- `sorted_indices`: Caller-provided permutation of `1:xnum` ordered by
    distance from `eruption_node`. Within a single `make_flow` call the
    eruption location is fixed, so the same permutation is reused across
    all pulses.
- `eruption_node`: Grid node index of the eruption. Caller is
    responsible for clamping into `1:xnum`.
- `eruption_thickness`: The thickness of lava at the eruption point (meters)
- `residual_lava_thickness`: The residual thickness of lava (meters)
- `tolerance`: The tolerance for the lava flow model
- `nmax`: The maximum number of iterations for the lava flow model

# Returns
- `icount`: The number of iterations for the lava flow model
"""
function lava_flow_pulse(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    lava_thickness::Vector{Float64},
    lava_thickness_old::Vector{Float64},
    sorted_indices::Vector{Int},
    eruption_node::Int,
    eruption_thickness::Float64,
    residual_lava_thickness::Float64;
    tolerance::Float64=1e-6,
    nmax::Int=20
)::Int
    xnum = length(topo_gridx)

    lava_thickness[eruption_node] += eruption_thickness
    fill!(lava_thickness_old, 1e38)

    max_diff = calculate_max_difference(lava_thickness, lava_thickness_old)
    icount = 0

    while max_diff > tolerance && icount < nmax
        lava_thickness_old .= lava_thickness
        @inbounds for mm in 1:xnum
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

                delta_elevation_left = max(0.0, y_surface_left - y_surface_active)
                delta_elevation_right = max(0.0, y_surface_right - y_surface_active)

                delta_elevation_nonzero_minimum = calculate_minimum_nonzero_delta_elevation(
                    delta_elevation_left, delta_elevation_right
                )

                delta_elevation_total = abs(delta_elevation_left) + abs(delta_elevation_right)

                if delta_elevation_total > 0.0 && thickness_active > residual_lava_thickness
                    total_potential_outflow_thickness = calculate_total_potential_outflow_thickness(
                        lava_thickness[i],
                        residual_lava_thickness,
                        delta_elevation_nonzero_minimum
                    )

                    total_outflow_thickness = 0.0
                    if delta_elevation_left > 0.0
                        outflow_thickness_to_left = calculate_outflow_thickness(
                            total_potential_outflow_thickness,
                            delta_elevation_left,
                            delta_elevation_total
                        )
                        total_outflow_thickness += outflow_thickness_to_left
                        lava_thickness[i-1] += outflow_thickness_to_left
                    end

                    if delta_elevation_right > 0.0
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

""" Radiate indices from a reference index, allocation-free.

Writes a permutation of `1:size` into `sorted_indices`, ordered by
ascending `abs(i - ref_index)`. Uses `distances_scratch` as a
preallocated scratch buffer. Bit-equivalent to `radiate_indices_old`
because `sortperm!` defaults to a stable sort: ties (same distance) are
broken by ascending original index, which is the same tie-breaker the
old `indices[sortperm(distances)]` formulation produces.

# Arguments
- `sorted_indices`: Output buffer, length `size`. Filled with the
    permutation; prior contents irrelevant.
- `distances_scratch`: Scratch buffer, length `size`. Overwritten.
- `size`: Length of the index range.
- `ref_index`: Reference index from which distance is measured.
"""
function radiate_indices!(
    sorted_indices::Vector{Int},
    distances_scratch::Vector{Int},
    size::Int,
    ref_index::Int
)::Nothing
    @inbounds for i in 1:size
        distances_scratch[i] = abs(i - ref_index)
    end
    sortperm!(sorted_indices, distances_scratch)
    return nothing
end

""" Original allocating implementation, kept for reference.

`collect(1:size)` + `abs.(... .- ref_index)` + `sortperm(distances)` +
`indices[sortperm(...)]` allocated four `Vector{Int}` of length `size`
per call; with `npulses ≈ 100` per flow these were a major component of
the lava-flow allocation budget. Replaced by `radiate_indices!` above.
"""
function radiate_indices_old(size::Int, ref_index::Int)::Vector{Int}
    indices = collect(1:size)
    distances = abs.(indices .- ref_index)
    return indices[sortperm(distances)]
end

""" Calculate the maximum difference between new and old lava thicknesses.

Allocation-free fused loop that replaces the original
`maximum(abs.(a .- b))` formulation. Bit-equivalent to that formulation
for finite inputs (the reduction is `max` over `abs` of element-wise
differences; `max` is associative so iteration order doesn't matter).
NaN propagation matches the original because `max(NaN, x) === NaN` and
`max(x, NaN) === NaN`.

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
    max_diff = 0.0
    @inbounds @simd for i in eachindex(lava_thickness)
        d = abs(lava_thickness[i] - lava_thickness_old[i])
        max_diff = max(max_diff, d)
    end
    return max_diff
end

""" Original allocating implementation, kept for reference.

The `abs.(...)` broadcast allocates a fresh `Vector{Float64}` of length
`xnum` every call, and `lava_flow_pulse` calls this function once per
iteration of its convergence loop (up to `nmax=1000` per pulse, with
many pulses per flow). Replaced by the fused-loop version above.
"""
function calculate_max_difference_old(
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