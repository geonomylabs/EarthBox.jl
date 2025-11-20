module GlobalIndices

""" Calculate global cell index for a given i and j.

# Arguments
- `i::Int`: The i-index of the cell (y-direction).
- `j::Int`: The j-index of the cell (x-direction).
- `ynum::Int`: The number basic grid nodes in the y-direction.

# Returns
- `cell_index::Int`: The global cell index.
"""
function get_global_basic_cell_index(i::Int, j::Int, ynum::Int)::Int
    return (j-1)*(ynum-1) + i
end

function get_global_ivx_unknown_index(cell_index::Int)::Int
    return cell_index*3 - 2
end

function get_global_ivy_unknown_index(ivx::Int)::Int
    return ivx + 1
end

function get_global_ipr_unknown_index(ivx::Int)::Int
    return ivx + 2
end

function get_cell_index_from_global_unknown(global_unknown_index::Int)::Int
    return floor(Int, (global_unknown_index + 2) / 3)
end

function get_ivx_from_global_unknown(global_unknown_index::Int)::Int
    return get_global_ivx_unknown_index(get_cell_index_from_global_unknown(global_unknown_index))
end

function get_ivy_from_global_unknown(global_unknown_index::Int)::Int
    return get_global_ivy_unknown_index(get_ivx_from_global_unknown(global_unknown_index))
end

function get_ipr_from_global_unknown(global_unknown_index::Int)::Int
    return get_global_ipr_unknown_index(get_ivx_from_global_unknown(global_unknown_index))
end

function get_ij_from_global_unknown(global_unknown_index::Int, ynum::Int)::Tuple{Int,Int}
    cell_index = get_cell_index_from_global_unknown(global_unknown_index)
    j = floor(Int, (cell_index - 1) / (ynum - 1)) + 1
    i = cell_index - (j - 1) * (ynum - 1)
    return i, j
end

end
