module InitializeMarkerRecycling

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs
using ..Sticky: get_sticky_material_ids
using ..Sticky: check_insticky

"""
    initialize_recycling(model::ModelData, inside_flags::Vector{Int8}) -> Tuple{Vector{Int64}, Int, Int}

Initialize recycling process for markers.

# Steps
1. Get indices of markers outside domain
2. Count markers to recycle
3. Print recycle count (if enabled)

# Arguments
- `model::ModelData`: Model data container

# Returns
- `outside_indices::Vector{Int64}`: Indices of markers outside the domain
- `nrecycle_air::Int`: Number of markers to recycle from sticky air/water
- `nrecycle_subsurface::Int`: Number of markers to recycle from lithosphere
"""
function initialize_recycling(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Tuple{Vector{Int64}, Int, Int}
    outside_indices = GridFuncs.get_indices_of_markers_outside_domain(model, inside_flags)
    nrecycle_air, nrecycle_subsurface = count_markers_to_recycle_opt(model, outside_indices)
    print_info = false
    if print_info
        print_recycle_count(nrecycle_air, nrecycle_subsurface)
    end
    return outside_indices, nrecycle_air, nrecycle_subsurface
end

"""
    count_markers_to_recycle(model::ModelData) -> Tuple{Int, Int}

Count the number of markers to recycle.

# Arguments
- `model::ModelData`: Model data object.

# Returns
- `nrecycle_air::Int`: Number of markers to recycle from sticky air/water
- `nrecycle_subsurface::Int`: Number of markers to recycle from lithosphere
"""
function count_markers_to_recycle(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Tuple{Int, Int}
    marker_matid = model.markers.arrays.material.marker_matid.array
    matids_sticky = get_sticky_material_ids(model)
    outside_indices = GridFuncs.get_indices_of_markers_outside_domain(model, inside_flags)
    noutside = length(outside_indices)
    nrecycle_air = 0
    nrecycle_subsurface = 0
    for ioutside in 1:noutside
        imarker = outside_indices[ioutside]
        matid = marker_matid[imarker]
        insticky = check_insticky(matid, matids_sticky)
        if insticky
            nrecycle_air += 1
        else
            nrecycle_subsurface += 1
        end
    end
    return nrecycle_air, nrecycle_subsurface
end

"""
    count_markers_to_recycle_opt(
        model::ModelData, outside_indices::Vector{Int64}) -> Tuple{Int, Int}

Count the number of markers to recycle.

# Arguments
- `model::ModelData`: Model data object.
- `outside_indices::Vector{Int64}`: Indices of markers outside the model domain.

# Returns
- `nrecycle_air::Int`: Number of markers to recycle from sticky air/water
- `nrecycle_subsurface::Int`: Number of markers to recycle from lithosphere
"""
function count_markers_to_recycle_opt(
    model::ModelData,
    outside_indices::Vector{Int64}
)::Tuple{Int, Int}
    marker_matid = model.markers.arrays.material.marker_matid.array
    matids_sticky = get_sticky_material_ids(model)
    noutside = length(outside_indices)
    nrecycle_air = 0
    nrecycle_subsurface = 0
    for ioutside in 1:noutside
        imarker = outside_indices[ioutside]
        matid = marker_matid[imarker]
        insticky = check_insticky(matid, matids_sticky)
        if insticky
            nrecycle_air += 1
        else
            nrecycle_subsurface += 1
        end
    end
    return nrecycle_air, nrecycle_subsurface
end

"""
    count_all_markers_to_recycle(model::ModelData) -> Int

Count the total number of markers to recycle regardless of type.

# Arguments
- `model::ModelData`: Model data object.

# Returns
- `Int`: Number of markers to recycle
"""
function count_all_markers_to_recycle(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Int
    outside_indices = GridFuncs.get_indices_of_markers_outside_domain(model, inside_flags)
    return length(outside_indices)
end

function print_recycle_count(
    nrecycle_air::Int,
    nrecycle_subsurface::Int
)
    println(
        ">> nrecycle_air: ", nrecycle_air,
        " nrecycle_subsurface: ", nrecycle_subsurface
    )
end


end # module 