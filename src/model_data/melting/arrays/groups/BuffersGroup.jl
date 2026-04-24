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
    Buffers(marknum::Int; nlayers::Int=50)

## Arguments
- `marknum::Int`: Number of markers (marker-sized buffers are sized to this length).
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

function Buffers(marknum::Int; nlayers::Int=50)::Buffers
    return Buffers(
        MarkerArrayInt1DState(
            zeros(Int64, marknum),
            ADATA.partial_melt_marker_indices.name,
            ADATA.partial_melt_marker_indices.units,
            ADATA.partial_melt_marker_indices.description
        ),
        MarkerArrayInt1DState(
            fill(Int64(-1), marknum),
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
            fill(Int64(-1), marknum),
            ADATA.layered_partial_melt_indices.name,
            ADATA.layered_partial_melt_indices.units,
            ADATA.layered_partial_melt_indices.description
        )
    )
end

end # module
