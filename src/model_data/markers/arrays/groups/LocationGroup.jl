"""
Module for marker location arrays.

Provides data structures for storing marker spatial coordinates.
"""
module LocationGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "location"

const ADATA = get_eb_arrays()

"""
    Location <: AbstractArrayGroup

Array group for marker spatial coordinates.

# Fields
- `marker_x::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_x.description)
- `marker_y::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_y.description)

# Nested Dot Access
- `marker_x = $(ROOT_NAME).$(GRP_NAME).marker_x.array`
- `marker_y = $(ROOT_NAME).$(GRP_NAME).marker_y.array`

# Constructor
    Location(marknum::Int)

Initializes marker location arrays with zero values.

## Arguments
- `marknum::Int`: Number of markers
"""
mutable struct Location <: AbstractArrayGroup
    marker_x::MarkerArrayFloat1DState{Float64}
    marker_y::MarkerArrayFloat1DState{Float64}
end

function Location(marknum::Int)::Location
    data = Location(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),     # array
            ADATA.marker_x.name,         # name
            ADATA.marker_x.units,        # units
            ADATA.marker_x.description   # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),     # array
            ADATA.marker_y.name,         # name
            ADATA.marker_y.units,        # units
            ADATA.marker_y.description   # description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(data::Location)
    data.marker_x.outform.header = "X_(m)"
    data.marker_y.outform.header = "Y_(m)"
    return nothing
end

end # module
