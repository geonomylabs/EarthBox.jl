"""
Module for marker solidification scratch buffers.

Holds pre-allocated marker-sized buffers used by
`MeltModel.Solidification.solidify!`. Currently a single field for per-
marker random numbers consumed by friction-coefficient randomization;
refilled via `Random.rand!` each call.
"""
module SolidificationGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "solidification"

const ADATA = get_eb_arrays()

"""
    Solidification <: AbstractArrayGroup

Array group holding scratch buffers reused by melt-solidification routines.

# Fields
- `marker_random_buffer::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_random_buffer.description)

# Nested Dot Access
- `model.markers.arrays.solidification.marker_random_buffer.array`

# Constructor
    Solidification(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers (buffer is sized to this length).
"""
mutable struct Solidification <: AbstractArrayGroup
    marker_random_buffer::MarkerArrayFloat1DState{Float64}
end

function Solidification(marknum::Int)::Solidification
    return Solidification(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.marker_random_buffer.name,
            ADATA.marker_random_buffer.units,
            ADATA.marker_random_buffer.description
        )
    )
end

end # module
