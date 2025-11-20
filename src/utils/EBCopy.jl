module EBCopy

""" Copy source array to target array.

# Arguments
- `source_array`: Source array
- `target_array`: Target array

# Updated Arrays
- `target_array`: Target array is updated with values from source array
"""
function copy_array_1d!(
    source_array::Vector{Float64},
    target_array::Vector{Float64}
)::Nothing
    if length(target_array) != length(source_array)
        error("Array sizes do not match.")
    end
    
    for i in 1:length(target_array)
        target_array[i] = source_array[i]
    end
    
    return nothing
end

""" Update elevation for vertical advection and diffusion.

# Arguments
- `topo_gridy_new`: New topography grid y-coordinates (meters)
- `gridt`: Multi-dimensional topography array. Y-coordinate (meters) of
    topography nodes are stored in gridt[2, :].

# Updated Arrays
- `gridt`: Multi-dimensional topography array is updated with new y-coordinates
"""
function copy_new_topography_to_topography_array(
    topo_gridy_new::Vector{Float64},
    gridt::Matrix{Float64}
)::Nothing
    toponum = length(topo_gridy_new)
    for i in 1:toponum
        gridt[2, i] = topo_gridy_new[i]
    end
    
    return nothing
end

function copy_topography_coordinate_arrays(gridt::Matrix{Float64})
    topo_gridx = copy(gridt[1, :])
    topo_gridy = copy(gridt[2, :])
    return topo_gridx, topo_gridy
end

end # module