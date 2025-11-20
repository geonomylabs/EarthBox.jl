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
function calculate_downstream_distances_for_nodes(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}}
    topo_shape = classify_topography_shape_at_each_node(topo_gridx, topo_gridy)
    drainage_divides_x = calculate_x_locations_of_drainage_divides(
        topo_gridx, topo_shape)
    toponum = length(topo_gridx)
    ndrainage = length(drainage_divides_x)
    downstream_distances = zeros(toponum)
    
    for i in 2:toponum
        xnode = topo_gridx[i]
        for j in 2:ndrainage
            x_left_drainage = drainage_divides_x[j-1]
            x_right_drainage = drainage_divides_x[j]
            if x_left_drainage <= xnode <= x_right_drainage
                distance = x_right_drainage - x_left_drainage
                downstream_distances[i] = distance
                break
            end
        end
    end
    return downstream_distances, drainage_divides_x
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
function classify_topography_shape_at_each_node(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64}
)::Vector{Int64}
    toponum = length(topo_gridx)
    topo_shape = fill(2, toponum)
    
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
    return topo_shape
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
function calculate_x_locations_of_drainage_divides(
    topo_gridx::Vector{Float64},
    topo_shape::Vector{Int64}
)::Vector{Float64}
    toponum = length(topo_gridx)
    drainage_divides_x_tmp = fill(-1e38, toponum)
    drainage_divides_x_tmp[1] = topo_gridx[1]
    icount = 1
    
    for i in 2:toponum
        if topo_shape[i] != 0
            drainage_divides_x_tmp[icount+1] = topo_gridx[i]
            icount += 1
        end
    end
    
    drainage_divides_x = zeros(icount)
    for i in 1:icount
        drainage_divides_x[i] = drainage_divides_x_tmp[i]
    end
    return drainage_divides_x
end

end # module 