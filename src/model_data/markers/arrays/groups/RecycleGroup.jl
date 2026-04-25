"""
Module for marker recycling scratch buffers.

Holds a pre-allocated marknum-sized buffer used by
`GridFuncs.get_indices_of_markers_outside_domain` to avoid allocating a
fresh `Vector{Int64}(undef, marknum)` on every recycling call.
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

# Nested Dot Access
- `model.markers.arrays.recycle.marker_outside_indices_scratch.array`

# Constructor
    Recycle(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers (buffer is sized to this length).
"""
mutable struct Recycle <: AbstractArrayGroup
    marker_outside_indices_scratch::MarkerArrayInt1DState{Int64}
end

function Recycle(marknum::Int)::Recycle
    return Recycle(
        MarkerArrayInt1DState(
            zeros(Int64, marknum),
            ADATA.marker_outside_indices_scratch.name,
            ADATA.marker_outside_indices_scratch.units,
            ADATA.marker_outside_indices_scratch.description
        )
    )
end

end # module
