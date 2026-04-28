module DownstreamDistance

""" Calculate downstream distances for each node.

# Arguments
- `topo_gridx`: Topography grid x-coordinates (meters)
- `topo_gridy`: Topography grid y-coordinates (meters). Note that y increases with
    depth

# Returns
- `downstream_distances`: Downstream distances for each node (meters). This 
    parameter represents the width of the drainage basin
- `drainage_divides_x`: X-locations of drainage divides (meters)
"""
function calculate_downstream_distances_for_nodes!(
    downstream_distances::Vector{Float64},
    drainage_divides_x::Vector{Float64},
    topo_shape_scratch::Vector{Int64},
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64}
)::Int
    classify_topography_shape_at_each_node!(topo_shape_scratch, topo_gridx, topo_gridy)
    ndivides = calculate_x_locations_of_drainage_divides!(
        drainage_divides_x, topo_gridx, topo_shape_scratch)
    toponum = length(topo_gridx)
    fill!(downstream_distances, 0.0)

    for i in 2:toponum
        xnode = topo_gridx[i]
        for j in 2:ndivides
            x_left_drainage = drainage_divides_x[j-1]
            x_right_drainage = drainage_divides_x[j]
            if x_left_drainage <= xnode <= x_right_drainage
                downstream_distances[i] = x_right_drainage - x_left_drainage
                break
            end
        end
    end
    return ndivides
end

""" Classify topography shape at each node.

# Arguments
- `topo_gridx`: Grid x-coordinates (meters)
- `topo_gridy`: Grid y-coordinates (meters). Note that y increases with depth

# Returns
- `topo_shape`: Flag for local minima and maxima:
    -1: local maximum node
     0: sloping node
     1: local minimum node
     2: flat node
"""
function classify_topography_shape_at_each_node!(
    topo_shape::Vector{Int64},
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64}
)::Nothing
    toponum = length(topo_gridx)
    fill!(topo_shape, 2)

    for i in 2:toponum-1
        y_left = topo_gridy[i-1]
        y_center = topo_gridy[i]
        y_right = topo_gridy[i+1]
        if y_center > y_left && y_center > y_right
            topo_shape[i] = 1
        elseif y_center < y_left && y_center < y_right
            topo_shape[i] = -1
        elseif y_center > y_right && y_center < y_left
            topo_shape[i] = 0
        elseif y_center < y_right && y_center > y_left
            topo_shape[i] = 0
        end
    end
    return nothing
end

""" Calculate x-locations of drainage divides.

# Arguments
- `topo_gridx`: Grid x-coordinates (meters)
- `topo_shape`: Flag for local minima and maxima:
    -1: local maximum node
     0: sloping node
     1: local minimum node
     2: flat node

# Returns
- `drainage_divides_x`: X-locations of drainage divides (meters)
"""
function calculate_x_locations_of_drainage_divides!(
    drainage_divides_x::Vector{Float64},
    topo_gridx::Vector{Float64},
    topo_shape::Vector{Int64}
)::Int
    toponum = length(topo_gridx)
    drainage_divides_x[1] = topo_gridx[1]
    icount = 1

    for i in 2:toponum
        if topo_shape[i] != 0
            drainage_divides_x[icount+1] = topo_gridx[i]
            icount += 1
        end
    end
    return icount
end

end # module 