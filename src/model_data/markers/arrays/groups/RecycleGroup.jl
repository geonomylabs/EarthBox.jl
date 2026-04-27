"""
Module for marker recycling scratch buffers.

Holds pre-allocated marknum-sized buffers used by
`GridFuncs.get_indices_of_markers_outside_domain` and
`GridFuncs.get_marker_inside_flags` to avoid per-call
`Vector{...}(undef, marknum)` allocations on every recycling/post-solver
call.
"""
module RecycleGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "recycle"

const ADATA = get_eb_arrays()

"""
    Recycle <: AbstractArrayGroup

Array group holding scratch buffers reused by marker-recycling routines.

# Fields
- `marker_outside_indices_scratch::`[`MarkerArrayInt1DState`](@ref) Int64: $(ADATA.marker_outside_indices_scratch.description)
- `marker_inside_flags_buffer::`[`MarkerArrayInt1DState`](@ref) Int8: $(ADATA.marker_inside_flags_buffer.description)

# Nested Dot Access
- `model.markers.arrays.recycle.marker_outside_indices_scratch.array`
- `model.markers.arrays.recycle.marker_inside_flags_buffer.array`

# Constructor
    Recycle(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers (buffers are sized to this length).
"""
mutable struct Recycle <: AbstractArrayGroup
    marker_outside_indices_scratch::MarkerArrayInt1DState{Int64}
    marker_inside_flags_buffer::MarkerArrayInt1DState{Int8}
end

function Recycle(marknum::Int)::Recycle
    return Recycle(
        MarkerArrayInt1DState(
            zeros(Int64, marknum),
            ADATA.marker_outside_indices_scratch.name,
            ADATA.marker_outside_indices_scratch.units,
            ADATA.marker_outside_indices_scratch.description
        ),
        MarkerArrayInt1DState(
            zeros(Int8, marknum),
            ADATA.marker_inside_flags_buffer.name,
            ADATA.marker_inside_flags_buffer.units,
            ADATA.marker_inside_flags_buffer.description
        )
    )
end

end # module
