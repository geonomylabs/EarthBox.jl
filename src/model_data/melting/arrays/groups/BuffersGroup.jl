module BuffersGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.melting.arrays"
const GRP_NAME = "buffers"

const ADATA = get_eb_arrays()

"""
    Buffers <: AbstractArrayGroup

Array group holding marker-sized scratch buffers reused across melt-extraction
calls, replacing per-call allocations. Each buffer follows the packed-prefix +
count convention: only positions `[1:count]` hold valid data after each call,
and the tail carries stale values from prior calls. Consumers must bound their
loops by the returned count, never by `length(buffer)`.

# Fields
- `partial_melt_marker_indices::`[`MarkerArrayInt1DState`](@ref) Int64: $(ADATA.partial_melt_marker_indices.description)
- `marker_indices_tmp::`[`MarkerArrayInt1DState`](@ref) Int64: $(ADATA.marker_indices_tmp.description)

# Nested Dot Access
- `partial_melt_marker_indices = $(ROOT_NAME).$(GRP_NAME).partial_melt_marker_indices.array`
- `marker_indices_tmp = $(ROOT_NAME).$(GRP_NAME).marker_indices_tmp.array`

# Constructor
    Buffers(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers (buffers are sized to this length).
"""
mutable struct Buffers <: AbstractArrayGroup
    partial_melt_marker_indices::MarkerArrayInt1DState{Int64}
    marker_indices_tmp::MarkerArrayInt1DState{Int64}
end

function Buffers(marknum::Int)::Buffers
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
        )
    )
end

end # module
