module BuffersGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.Arrays.ArrayTypes.Array1DInt: Array1DIntState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.melting.arrays"
const GRP_NAME = "buffers"

const ADATA = get_eb_arrays()

"""
    Buffers <: AbstractArrayGroup

Array group holding scratch buffers reused across melt-extraction calls,
replacing per-call allocations. Each buffer follows the packed-prefix +
count convention: only positions `[1:count]` hold valid data after each call,
and the tail carries stale values from prior calls. Consumers must bound their
loops by the returned count, never by `length(buffer)`.

The three marker-sized fields (`partial_melt_marker_indices`,
`marker_indices_tmp`, `layered_partial_melt_indices`) are only consumed inside
the melt-extraction code path (gated by `iuse_melting == 1 && iuse_extraction
== 1`). They are constructed at length 0 and lazily resized to `marknum` on
the first call to `update_melt_extraction!` via
[`ensure_extraction_buffers!`](@ref) — keeping RAM low when extraction is off.

# Fields
- `partial_melt_marker_indices::`[`MarkerArrayInt1DState`](@ref) Int64: $(ADATA.partial_melt_marker_indices.description)
- `marker_indices_tmp::`[`MarkerArrayInt1DState`](@ref) Int64: $(ADATA.marker_indices_tmp.description)
- `pm_layer_counts::`[`Array1DIntState`](@ref): $(ADATA.pm_layer_counts.description)
- `pm_layer_offsets::`[`Array1DIntState`](@ref): $(ADATA.pm_layer_offsets.description)
- `layered_partial_melt_indices::`[`MarkerArrayInt1DState`](@ref) Int64: $(ADATA.layered_partial_melt_indices.description)

# Nested Dot Access
- `partial_melt_marker_indices = $(ROOT_NAME).$(GRP_NAME).partial_melt_marker_indices.array`
- `marker_indices_tmp = $(ROOT_NAME).$(GRP_NAME).marker_indices_tmp.array`
- `pm_layer_counts = $(ROOT_NAME).$(GRP_NAME).pm_layer_counts.array`
- `pm_layer_offsets = $(ROOT_NAME).$(GRP_NAME).pm_layer_offsets.array`
- `layered_partial_melt_indices = $(ROOT_NAME).$(GRP_NAME).layered_partial_melt_indices.array`

# Constructor
    Buffers(; nlayers::Int=50)

The marker-length buffers start empty and are sized on first use by
[`ensure_extraction_buffers!`](@ref).

## Arguments
- `nlayers::Int`: Number of layers used to bin partially molten mantle markers
  for layered shallowest-marker search. Must match the nlayers used in
  `construct_layered_partially_molten_arrays!`. Derived at runtime from
  `length(pm_layer_counts)`.
"""
mutable struct Buffers <: AbstractArrayGroup
    partial_melt_marker_indices::MarkerArrayInt1DState{Int64}
    marker_indices_tmp::MarkerArrayInt1DState{Int64}
    pm_layer_counts::Array1DIntState
    pm_layer_offsets::Array1DIntState
    layered_partial_melt_indices::MarkerArrayInt1DState{Int64}
end

function Buffers(; nlayers::Int=50)::Buffers
    return Buffers(
        MarkerArrayInt1DState(
            Int64[],
            ADATA.partial_melt_marker_indices.name,
            ADATA.partial_melt_marker_indices.units,
            ADATA.partial_melt_marker_indices.description
        ),
        MarkerArrayInt1DState(
            Int64[],
            ADATA.marker_indices_tmp.name,
            ADATA.marker_indices_tmp.units,
            ADATA.marker_indices_tmp.description
        ),
        Array1DIntState(
            zeros(Int, nlayers),
            ADATA.pm_layer_counts.name,
            ADATA.pm_layer_counts.units,
            ADATA.pm_layer_counts.description
        ),
        Array1DIntState(
            zeros(Int, nlayers + 1),
            ADATA.pm_layer_offsets.name,
            ADATA.pm_layer_offsets.units,
            ADATA.pm_layer_offsets.description
        ),
        MarkerArrayInt1DState(
            Int64[],
            ADATA.layered_partial_melt_indices.name,
            ADATA.layered_partial_melt_indices.units,
            ADATA.layered_partial_melt_indices.description
        )
    )
end

"""
    ensure_extraction_buffers!(buffers::Buffers, marknum::Int)

Idempotently size the marker-length extraction scratch buffers to `marknum`
and reproduce the original constructor's initial values on the transition
from length 0. Safe to call every timestep — does nothing when sizes already
match. Called at the top of `MeltModel.Extraction.update_melt_extraction!`
when `iuse_melting == 1 && iuse_extraction == 1`.
"""
function ensure_extraction_buffers!(buffers::Buffers, marknum::Int)::Nothing
    _ensure_buffer!(buffers.partial_melt_marker_indices.array, marknum, Int64(0))
    _ensure_buffer!(buffers.marker_indices_tmp.array, marknum, Int64(-1))
    _ensure_buffer!(buffers.layered_partial_melt_indices.array, marknum, Int64(-1))
    return nothing
end

function _ensure_buffer!(buf::Vector{Int64}, target::Int, sentinel::Int64)::Nothing
    current = length(buf)
    if current == target
        return nothing
    end
    resize!(buf, target)
    if target > current
        @inbounds for i in (current + 1):target
            buf[i] = sentinel
        end
    end
    return nothing
end

end # module
