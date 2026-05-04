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

Each buffer is constructed at length 0 and lazily resized to `marknum`
on first use via [`ensure_compaction_buffers!`](@ref). All nine fields
are flagged with `ibackup=false` so the JLD2 backup loader does not try
to restore them — they are scratch and are recomputed on each call.
The lazy-init hook is invoked at the top of:
- `ApplyCompaction.apply_compaction_model!` (lava flow, salt deposition)
- The `compaction_correction_type == "variable_property"` branch in
  `SedimentTransportSolverManager.run_sediment_transport_time_steps!`
  (downhill-diffusion sediment transport with compaction correction)
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
    Compaction()

The marker-length buffers start empty and are sized on first use by
[`ensure_compaction_buffers!`](@ref).
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

function Compaction()::Compaction
    return Compaction(
        MarkerArrayInt1DState(
            Int64[],
            ADATA.markers_topo_xindex_buffer.name,
            ADATA.markers_topo_xindex_buffer.units,
            ADATA.markers_topo_xindex_buffer.description;
            ibackup=false
        ),
        MarkerArrayInt1DState(
            Int64[],
            ADATA.markers_compaction_yindex_buffer.name,
            ADATA.markers_compaction_yindex_buffer.units,
            ADATA.markers_compaction_yindex_buffer.description;
            ibackup=false
        ),
        MarkerArrayFloat1DState(
            Float64[],
            ADATA.markers_unit_distance_from_cell_top_buffer.name,
            ADATA.markers_unit_distance_from_cell_top_buffer.units,
            ADATA.markers_unit_distance_from_cell_top_buffer.description;
            ibackup=false
        ),
        MarkerArrayFloat1DState(
            Float64[],
            ADATA.total_marker_compaction_displacement_buffer.name,
            ADATA.total_marker_compaction_displacement_buffer.units,
            ADATA.total_marker_compaction_displacement_buffer.description;
            ibackup=false
        ),
        MarkerArrayFloat1DState(
            Float64[],
            ADATA.sticky_displacement_factors_buffer.name,
            ADATA.sticky_displacement_factors_buffer.units,
            ADATA.sticky_displacement_factors_buffer.description;
            ibackup=false
        ),
        MarkerArrayFloat1DState(
            Float64[],
            ADATA.sticky_marker_displacement_buffer.name,
            ADATA.sticky_marker_displacement_buffer.units,
            ADATA.sticky_marker_displacement_buffer.description;
            ibackup=false
        ),
        MarkerArrayInt1DState(
            Int64[],
            ADATA.marker_swarm_index_scratch.name,
            ADATA.marker_swarm_index_scratch.units,
            ADATA.marker_swarm_index_scratch.description;
            ibackup=false
        ),
        MarkerArrayFloat1DState(
            Float64[],
            ADATA.marker_swarm_x_gather_scratch.name,
            ADATA.marker_swarm_x_gather_scratch.units,
            ADATA.marker_swarm_x_gather_scratch.description;
            ibackup=false
        ),
        MarkerArrayInt1DState(
            Int64[],
            ADATA.marker_swarm_sortperm_scratch.name,
            ADATA.marker_swarm_sortperm_scratch.units,
            ADATA.marker_swarm_sortperm_scratch.description;
            ibackup=false
        ),
    )
end

"""
    ensure_compaction_buffers!(compaction::Compaction, marknum::Int)

Idempotently size all nine marker-length compaction scratch buffers to
`marknum`. Safe to call repeatedly — does nothing when sizes already
match. Called at the top of:
- `ApplyCompaction.apply_compaction_model!`
- The `variable_property` branch in
  `SedimentTransportSolverManager.run_sediment_transport_time_steps!`
"""
function ensure_compaction_buffers!(
    compaction::Compaction, marknum::Int
)::Nothing
    _ensure_int_buffer!(compaction.markers_topo_xindex_buffer.array, marknum)
    _ensure_int_buffer!(compaction.markers_compaction_yindex_buffer.array, marknum)
    _ensure_float_buffer!(compaction.markers_unit_distance_from_cell_top_buffer.array, marknum)
    _ensure_float_buffer!(compaction.total_marker_compaction_displacement_buffer.array, marknum)
    _ensure_float_buffer!(compaction.sticky_displacement_factors_buffer.array, marknum)
    _ensure_float_buffer!(compaction.sticky_marker_displacement_buffer.array, marknum)
    _ensure_int_buffer!(compaction.marker_swarm_index_scratch.array, marknum)
    _ensure_float_buffer!(compaction.marker_swarm_x_gather_scratch.array, marknum)
    _ensure_int_buffer!(compaction.marker_swarm_sortperm_scratch.array, marknum)
    return nothing
end

function _ensure_int_buffer!(buf::Vector{Int64}, target::Int)::Nothing
    current = length(buf)
    if current == target
        return nothing
    end
    resize!(buf, target)
    if target > current
        @inbounds for i in (current + 1):target
            buf[i] = Int64(0)
        end
    end
    return nothing
end

function _ensure_float_buffer!(buf::Vector{Float64}, target::Int)::Nothing
    current = length(buf)
    if current == target
        return nothing
    end
    resize!(buf, target)
    if target > current
        @inbounds for i in (current + 1):target
            buf[i] = 0.0
        end
    end
    return nothing
end

end # module
