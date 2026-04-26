"""
Module for marker compaction scratch buffers.

Holds pre-allocated marker-sized buffers used internally by
`MarkerCompaction.compact_sediment_and_advect_markers` and its helpers
(sticky path included). All buffers are single-use within one compaction
call and never shared with any other call site or cross-domain consumer.

The swarm-index / sortperm buffers below are also reused across the
back-to-back sedimentary-basin and sticky passes inside
`calculate_swarm_indices_for_sediment_and_sticky` and
`calculate_x_sorted_swarm_indices`; both passes are sequential on the
main thread and each writes its result into a fresh tight `Vector{Int64}`
before the next pass starts, so a single shared scratch per role is safe.

Each buffer is sized to `marknum` at construction. Production callers
pass these via positional arguments to `compact_sediment_and_advect_markers`
so the function signature stays explicit (no `Union{T, Nothing}=nothing`
keyword-arg dispatch hazard).
"""
module CompactionGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "compaction"

const ADATA = get_eb_arrays()

"""
    Compaction <: AbstractArrayGroup

Array group holding scratch buffers used by marker compaction.

# Fields
- `markers_topo_xindex_buffer::`[`MarkerArrayInt1DState`](@ref) Int64
- `markers_compaction_yindex_buffer::`[`MarkerArrayInt1DState`](@ref) Int64
- `markers_unit_distance_from_cell_top_buffer::`[`MarkerArrayFloat1DState`](@ref) Float64
- `total_marker_compaction_displacement_buffer::`[`MarkerArrayFloat1DState`](@ref) Float64
- `sticky_displacement_factors_buffer::`[`MarkerArrayFloat1DState`](@ref) Float64
- `sticky_marker_displacement_buffer::`[`MarkerArrayFloat1DState`](@ref) Float64
- `marker_swarm_index_scratch::`[`MarkerArrayInt1DState`](@ref) Int64
- `marker_swarm_x_gather_scratch::`[`MarkerArrayFloat1DState`](@ref) Float64
- `marker_swarm_sortperm_scratch::`[`MarkerArrayInt1DState`](@ref) Int64

# Constructor
    Compaction(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers (each buffer is sized to this length).
"""
mutable struct Compaction <: AbstractArrayGroup
    markers_topo_xindex_buffer::MarkerArrayInt1DState{Int64}
    markers_compaction_yindex_buffer::MarkerArrayInt1DState{Int64}
    markers_unit_distance_from_cell_top_buffer::MarkerArrayFloat1DState{Float64}
    total_marker_compaction_displacement_buffer::MarkerArrayFloat1DState{Float64}
    sticky_displacement_factors_buffer::MarkerArrayFloat1DState{Float64}
    sticky_marker_displacement_buffer::MarkerArrayFloat1DState{Float64}
    marker_swarm_index_scratch::MarkerArrayInt1DState{Int64}
    marker_swarm_x_gather_scratch::MarkerArrayFloat1DState{Float64}
    marker_swarm_sortperm_scratch::MarkerArrayInt1DState{Int64}
end

function Compaction(marknum::Int)::Compaction
    return Compaction(
        MarkerArrayInt1DState(
            zeros(Int64, marknum),
            ADATA.markers_topo_xindex_buffer.name,
            ADATA.markers_topo_xindex_buffer.units,
            ADATA.markers_topo_xindex_buffer.description
        ),
        MarkerArrayInt1DState(
            zeros(Int64, marknum),
            ADATA.markers_compaction_yindex_buffer.name,
            ADATA.markers_compaction_yindex_buffer.units,
            ADATA.markers_compaction_yindex_buffer.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.markers_unit_distance_from_cell_top_buffer.name,
            ADATA.markers_unit_distance_from_cell_top_buffer.units,
            ADATA.markers_unit_distance_from_cell_top_buffer.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.total_marker_compaction_displacement_buffer.name,
            ADATA.total_marker_compaction_displacement_buffer.units,
            ADATA.total_marker_compaction_displacement_buffer.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.sticky_displacement_factors_buffer.name,
            ADATA.sticky_displacement_factors_buffer.units,
            ADATA.sticky_displacement_factors_buffer.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.sticky_marker_displacement_buffer.name,
            ADATA.sticky_marker_displacement_buffer.units,
            ADATA.sticky_marker_displacement_buffer.description
        ),
        MarkerArrayInt1DState(
            zeros(Int64, marknum),
            ADATA.marker_swarm_index_scratch.name,
            ADATA.marker_swarm_index_scratch.units,
            ADATA.marker_swarm_index_scratch.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.marker_swarm_x_gather_scratch.name,
            ADATA.marker_swarm_x_gather_scratch.units,
            ADATA.marker_swarm_x_gather_scratch.description
        ),
        MarkerArrayInt1DState(
            zeros(Int64, marknum),
            ADATA.marker_swarm_sortperm_scratch.name,
            ADATA.marker_swarm_sortperm_scratch.units,
            ADATA.marker_swarm_sortperm_scratch.description
        ),
    )
end

end # module
