module PartiallyMoltenZone

""" Populate pre-allocated layer buffers with partially molten marker indices.

Bins partially molten mantle markers (whose material id is in
`mantle_emplacement_mat_ids`) into `nlayers` equal y-range slabs by counting
sort, writing:
- `layer_counts[i]`: number of valid markers in layer `i`
- `layer_offsets[i]`: exclusive prefix-sum so layer `i`'s indices occupy
  `layered_partial_melt_indices[layer_offsets[i]+1 : layer_offsets[i+1]]`
- `layered_partial_melt_indices[layer_offsets[i]+1 : layer_offsets[i+1]]`:
  packed marker indices belonging to layer `i`

The number of layers is taken from `length(layer_counts)`. The function
performs three O(nmarkers_partial_melt) passes (y-range scan, count, pack)
and zero heap allocations.

`layer_offsets` is the exclusive prefix sum of `layer_counts`: `layer_offsets[1]
= 0` and `layer_offsets[i+1] = layer_offsets[i] + layer_counts[i]`. It locates
each layer's block inside the flat buffer, so layer `i`'s valid marker indices
are `layered_partial_melt_indices[layer_offsets[i]+1 : layer_offsets[i+1]]`
(count `layer_counts[i]`). This replaces a `Vector{Vector{Int64}}` with one
contiguous vector plus bookkeeping.

# Arguments
- `marker_x::Vector{Float64}`: Marker x-coordinates.
- `marker_y::Vector{Float64}`: Marker y-coordinates.
- `marker_matid::Vector{Int16}`: Marker material ids.
- `nmarkers_partial_melt::Int`: Number of valid entries in
  `partial_melt_marker_indices`.
- `mantle_emplacement_mat_ids::Vector{Int16}`: Material ids considered
  partially molten mantle for this call.
- `partial_melt_marker_indices::Vector{Int64}`: Candidate marker indices
  (only the first `nmarkers_partial_melt` entries are read).
- `layer_counts::Vector{Int}`: Scratch buffer of length `nlayers`; populated
  with per-layer counts.
- `layer_offsets::Vector{Int}`: Scratch buffer of length `nlayers + 1`;
  populated with exclusive prefix-sum offsets.
- `layered_partial_melt_indices::Vector{Int64}`: Scratch buffer packed with
  marker indices grouped by layer. Must be large enough to hold every marker
  that passes the matid filter (typically sized to `marknum`).
"""
function construct_layered_partially_molten_arrays!(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    nmarkers_partial_melt::Int,
    mantle_emplacement_mat_ids::Vector{Int16},
    partial_melt_marker_indices::Vector{Int64},
    layer_counts::Vector{Int},
    layer_offsets::Vector{Int},
    layered_partial_melt_indices::Vector{Int64}
)::Nothing
    nlayers = length(layer_counts)

    _, _, pm_ymin, pm_ymax = get_limits_of_partially_molten_zone(
        marker_x, marker_y, marker_matid, nmarkers_partial_melt,
        mantle_emplacement_mat_ids, partial_melt_marker_indices
    )

    # dy_inv = nlayers / (pm_ymax - pm_ymin); zero when the zone is degenerate
    # so every valid marker falls into layer 1 (matches old behavior where
    # ymin == ymax made all layers coincide).
    has_range = pm_ymax > pm_ymin
    dy_inv = has_range ? nlayers / (pm_ymax - pm_ymin) : 0.0

    fill!(layer_counts, 0)
    for j in 1:nmarkers_partial_melt
        imarker = partial_melt_marker_indices[j]
        matid = marker_matid[imarker]
        if matid in mantle_emplacement_mat_ids
            if has_range
                y = marker_y[imarker]
                ilayer = floor(Int, (y - pm_ymin) * dy_inv) + 1
                if ilayer < 1
                    ilayer = 1
                elseif ilayer > nlayers
                    ilayer = nlayers
                end
            else
                ilayer = 1
            end
            layer_counts[ilayer] += 1
        end
    end

    layer_offsets[1] = 0
    for i in 1:nlayers
        layer_offsets[i + 1] = layer_offsets[i] + layer_counts[i]
    end

    fill!(layer_counts, 0)
    for j in 1:nmarkers_partial_melt
        imarker = partial_melt_marker_indices[j]
        matid = marker_matid[imarker]
        if matid in mantle_emplacement_mat_ids
            if has_range
                y = marker_y[imarker]
                ilayer = floor(Int, (y - pm_ymin) * dy_inv) + 1
                if ilayer < 1
                    ilayer = 1
                elseif ilayer > nlayers
                    ilayer = nlayers
                end
            else
                ilayer = 1
            end
            pos = layer_offsets[ilayer] + layer_counts[ilayer] + 1
            layered_partial_melt_indices[pos] = imarker
            layer_counts[ilayer] += 1
        end
    end

    return nothing
end

""" Get limits of partially molten zone.

# Returns
- xmin::Float64
    - Minimum x-coordinate of the partially molten zone.
- xmax::Float64
    - Maximum x-coordinate of the partially molten zone.
- ymin::Float64
    - Minimum y-coordinate of the partially molten zone.
- ymax::Float64
    - Maximum y-coordinate of the partially molten zone.
"""
function get_limits_of_partially_molten_zone(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    nmarkers_partial_melt::Int,
    mantle_emplacement_mat_ids::Vector{Int16},
    partial_melt_marker_indices::Vector{Int64}
)::Tuple{Float64, Float64, Float64, Float64}
    xmin = 1e39
    xmax = -1e39
    ymin = 1e39
    ymax = -1e39
    # Get dimensions of partially molten zone
    for j in 1:nmarkers_partial_melt
        # get index of partially molten mantle marker
        imarker = partial_melt_marker_indices[j]
        y_marker = marker_y[imarker]
        x_marker = marker_x[imarker]
        matid = marker_matid[imarker]
        # Only consider ID's associated with the mantle melting model
        if matid in mantle_emplacement_mat_ids
            if y_marker < ymin
                ymin = y_marker
            end
            if y_marker > ymax
                ymax = y_marker
            end
            if x_marker < xmin
                xmin = x_marker
            end
            if x_marker > xmax
                xmax = x_marker
            end
        end
    end
    debug = false
    if debug
        println("   >> Partially molten zone limits: ", xmin, " ", xmax, " ", ymin, " ", ymax)
    end
    return xmin, xmax, ymin, ymax
end

end # module
