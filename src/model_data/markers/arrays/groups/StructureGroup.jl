"""
Module for marker structure-finding scratch buffers.

Holds pre-allocated marker-sized buffers used by layer-finding utilities such
as `TopAndBottom.filter_markers_outside_of_layer!`, called from drainage,
fractionation, and topography construction. Each buffer follows the
packed-prefix + count convention: only positions `[1:count]` hold valid data
after each call. Consumers must bound their loops by the returned count, never
by `length(buffer)`.
"""
module StructureGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "structure"

const ADATA = get_eb_arrays()

"""
    Structure <: AbstractArrayGroup

Array group holding scratch buffers reused across calls to layer-finding
utilities. Replaces per-call `Vector{Int64}(undef, marknum)` allocations in
`filter_markers_outside_of_layer`.

# Fields
- `marker_indices_layer::`[`MarkerArrayInt1DState`](@ref) Int64: $(ADATA.marker_indices_layer.description)

# Nested Dot Access
- `marker_indices_layer = $(ROOT_NAME).$(GRP_NAME).marker_indices_layer.array`

# Constructor
    Structure(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers (buffer is sized to this length).
"""
mutable struct Structure <: AbstractArrayGroup
    marker_indices_layer::MarkerArrayInt1DState{Int64}
end

function Structure(marknum::Int)::Structure
    return Structure(
        MarkerArrayInt1DState(
            zeros(Int64, marknum),
            ADATA.marker_indices_layer.name,
            ADATA.marker_indices_layer.units,
            ADATA.marker_indices_layer.description
        )
    )
end

end # module
