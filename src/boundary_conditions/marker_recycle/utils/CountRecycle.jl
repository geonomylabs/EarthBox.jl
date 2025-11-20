module CountRecycle

using EarthBox.ModelDataContainer: ModelData
using EarthBox: GridFuncs

"""
    get_left_and_right_recycle_counts(nrecycle::Int; use_asymmetric_extension::Bool=false) -> Tuple{Int, Int}

Get left and right recycle counts for extension models.

# Arguments
- `nrecycle::Int`: Total number of markers to recycle
- `use_asymmetric_extension::Bool=false`: Whether to use asymmetric extension

# Returns
- `Tuple{Int, Int}`: Number of markers to recycle on left and right sides
"""
function get_left_and_right_recycle_counts(
    nrecycle::Int;
    use_asymmetric_extension::Bool=false
)::Tuple{Int, Int}
    if use_asymmetric_extension
        nrecycle_left = 0
    else
        nrecycle_left = floor(Int, nrecycle/2)
    end
    nrecycle_right = nrecycle - nrecycle_left
    return nrecycle_left, nrecycle_right
end

end # module