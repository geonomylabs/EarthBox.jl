module PartiallyMoltenZone

""" Calculate partially molten marker indices and counts for each layer.

# Returns
- layer_counts::Vector{Int64}
    - Number of partially molten markers in each layer.
- layered_partially_molten_marker_indices::Vector{Vector{Int64}}
    - Indices of partially molten markers in each layer.
"""
function construct_layered_partially_molten_arrays(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    nmarkers_partial_melt::Int,
    mantle_emplacement_mat_ids::Vector{Int16},
    partial_melt_marker_indices::Vector{Int64};
    nlayers::Int=50
)::Tuple{Vector{Int64}, Vector{Vector{Int64}}}
    pm_xmin, pm_xmax, pm_ymin, pm_ymax = get_limits_of_partially_molten_zone(
        marker_x, marker_y, marker_matid, nmarkers_partial_melt,
        mantle_emplacement_mat_ids, partial_melt_marker_indices
    )

    layer_dims = create_dimension_array_for_layers(
        nlayers, pm_xmin, pm_xmax, pm_ymin, pm_ymax
    )

    layer_counts, layered_partially_molten_marker_indices = 
        calculate_partially_molten_marker_indices_for_layers(
            marker_x, marker_y, marker_matid, nmarkers_partial_melt,
            mantle_emplacement_mat_ids, partial_melt_marker_indices, layer_dims
        )

    return layer_counts, layered_partially_molten_marker_indices
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

""" Create array to store dimensions of each layers.

For each layer store: xmin, xmax, ymin, ymax.

# Returns
- layer_dims::Matrix{Float64}
    - Array to store dimensions of each layer:
    - layer_dims[i, 1] = xmin
    - layer_dims[i, 2] = xmax
    - layer_dims[i, 3] = ymin
    - layer_dims[i, 4] = ymax
"""
function create_dimension_array_for_layers(
    nlayers::Int,
    pm_xmin::Float64,
    pm_xmax::Float64,
    pm_ymin::Float64,
    pm_ymax::Float64
)::Matrix{Float64}
    dy = (pm_ymax - pm_ymin) / nlayers
    layer_dims = Matrix{Float64}(undef, nlayers, 4) #zeros(Float64, nlayers, 4)
    for i in 1:nlayers
        layer_dims[i, 1] = pm_xmin
        layer_dims[i, 2] = pm_xmax

        ytop = pm_ymin + (i-1)*dy
        layer_dims[i, 3] = ytop
        layer_dims[i, 4] = ytop + dy
    end
    return layer_dims
end

""" Calculate indices of partially molten markers in layers.

A multi-dimensional array is calculated to store the indices of partially
molten markers in layers created by dividing up the partially molten mantle 
domain into nlayers with dimensions stored in layer_dims. The number of 
markers in each layer is also calculated and stored in a separate array 
called layer_counts.

Only markers with that are partially molten and have material IDs that are
in mantle_emplacement_mat_ids are considered.

# Returns
- layer_counts::Vector{Int64}
    - Number of partially molten markers in each layer. This number is equal 
    to or less than nmarkers_partial_melt.
- layered_partially_molten_marker_indices::Vector{Vector{Int64}}
    - Array to store IDs of partially molten markers in each layer.
"""
function calculate_partially_molten_marker_indices_for_layers(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    marker_matid::Vector{Int16},
    nmarkers_partial_melt::Int,
    mantle_emplacement_mat_ids::Vector{Int16},
    partial_melt_flags::Vector{Int64},
    layer_dims::Matrix{Float64}
)::Tuple{Vector{Int64}, Vector{Vector{Int64}}}
    nlayers = size(layer_dims, 1)
    layer_counts = Vector{Int64}(undef, nlayers)

    layered_partially_molten_marker_indices = 
        create_layer_array_for_partially_molten_ids(nlayers, nmarkers_partial_melt)

    for i in 1:nlayers
        xmin = layer_dims[i, 1]
        xmax = layer_dims[i, 2]
        ymin = layer_dims[i, 3]
        ymax = layer_dims[i, 4]
        icount = 0
        for j in 1:nmarkers_partial_melt
            imarker = partial_melt_flags[j]
            y_marker = marker_y[imarker]
            x_marker = marker_x[imarker]
            matid = marker_matid[imarker]
            if matid in mantle_emplacement_mat_ids
                if xmin <= x_marker <= xmax && ymin <= y_marker <= ymax
                    layered_partially_molten_marker_indices[i][icount + 1] = imarker
                    icount += 1
                end
            end
        end
        layer_counts[i] = icount
    end
    return layer_counts, layered_partially_molten_marker_indices
end

""" Create array to store indices of partially molten markers in each layer.

# Returns
- layered_partially_molten_marker_indices::Vector{Vector{Int64}}
    - Indices of partially molten markers in each layer.
"""
function create_layer_array_for_partially_molten_ids(
    nlayers::Int,
    nmarkers_partial_melt::Int
)::Vector{Vector{Int64}}
    return [fill(-1, nmarkers_partial_melt) for _ in 1:nlayers]
end

end # module 