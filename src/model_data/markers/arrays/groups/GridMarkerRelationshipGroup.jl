"""
Module for grid-marker relationship arrays.

Provides data structures for storing the spatial relationship between markers
and grid nodes, including indices and normalized distances for interpolation.
"""
module GridMarkerRelationshipGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "grid_marker_relationship"

const ADATA = get_eb_arrays()

"""
    GridMarkerRelationship <: AbstractArrayGroup

Array group for grid-marker spatial relationships.

# Fields
- `marker_xn::`[`MarkerArrayInt1DState`](@ref) Int32: $(ADATA.marker_xn.description)
- `marker_yn::`[`MarkerArrayInt1DState`](@ref) Int32: $(ADATA.marker_yn.description)
- `marker_dx::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_dx.description)
- `marker_dy::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_dy.description)
- `marker_xn_vy::`[`MarkerArrayInt1DState`](@ref) Int32: $(ADATA.marker_xn_vy.description)
- `marker_yn_vy::`[`MarkerArrayInt1DState`](@ref) Int32: $(ADATA.marker_yn_vy.description)
- `marker_dx_vy::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_dx_vy.description)
- `marker_dy_vy::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_dy_vy.description)

# Nested Dot Access
- `marker_xn = $(ROOT_NAME).$(GRP_NAME).marker_xn.array`
- `marker_yn = $(ROOT_NAME).$(GRP_NAME).marker_yn.array`
- `marker_dx = $(ROOT_NAME).$(GRP_NAME).marker_dx.array`
- `marker_dy = $(ROOT_NAME).$(GRP_NAME).marker_dy.array`
- `marker_xn_vy = $(ROOT_NAME).$(GRP_NAME).marker_xn_vy.array`
- `marker_yn_vy = $(ROOT_NAME).$(GRP_NAME).marker_yn_vy.array`
- `marker_dx_vy = $(ROOT_NAME).$(GRP_NAME).marker_dx_vy.array`
- `marker_dy_vy = $(ROOT_NAME).$(GRP_NAME).marker_dy_vy.array`

# Constructor
    GridMarkerRelationship(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers
"""
mutable struct GridMarkerRelationship <: AbstractArrayGroup
    marker_xn::MarkerArrayInt1DState{Int32}
    marker_yn::MarkerArrayInt1DState{Int32}
    marker_dx::MarkerArrayFloat1DState{Float64}
    marker_dy::MarkerArrayFloat1DState{Float64}
    marker_xn_vy::MarkerArrayInt1DState{Int32}
    marker_yn_vy::MarkerArrayInt1DState{Int32}
    marker_dx_vy::MarkerArrayFloat1DState{Float64}
    marker_dy_vy::MarkerArrayFloat1DState{Float64}
end

function GridMarkerRelationship(marknum::Int)::GridMarkerRelationship
    data = GridMarkerRelationship(
        MarkerArrayInt1DState(
            zeros(Int32, marknum),           # array
            ADATA.marker_xn.name,            # name
            ADATA.marker_xn.units,           # units
            ADATA.marker_xn.description      # description
        ),
        MarkerArrayInt1DState(
            zeros(Int32, marknum),           # array
            ADATA.marker_yn.name,            # name
            ADATA.marker_yn.units,           # units
            ADATA.marker_yn.description      # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),         # array
            ADATA.marker_dx.name,            # name
            ADATA.marker_dx.units,           # units
            ADATA.marker_dx.description      # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),         # array
            ADATA.marker_dy.name,            # name
            ADATA.marker_dy.units,           # units
            ADATA.marker_dy.description      # description
        ),
        MarkerArrayInt1DState(
            zeros(Int32, marknum),           # array
            ADATA.marker_xn_vy.name,         # name
            ADATA.marker_xn_vy.units,        # units
            ADATA.marker_xn_vy.description   # description
        ),
        MarkerArrayInt1DState(
            zeros(Int32, marknum),           # array
            ADATA.marker_yn_vy.name,         # name
            ADATA.marker_yn_vy.units,        # units
            ADATA.marker_yn_vy.description   # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),         # array
            ADATA.marker_dx_vy.name,         # name
            ADATA.marker_dx_vy.units,        # units
            ADATA.marker_dx_vy.description   # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),         # array
            ADATA.marker_dy_vy.name,         # name
            ADATA.marker_dy_vy.units,        # units
            ADATA.marker_dy_vy.description   # description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(
    data::GridMarkerRelationship
)::Nothing
    data.marker_xn.outform.header = "UpperLeftBasicNodeIndex_xn"
    data.marker_yn.outform.header = "UpperLeftBasicNodeIndex_yn"
    data.marker_dx.outform.header = "dx_upperleft_basic"
    data.marker_dy.outform.header = "dy_upperleft_basic"
    data.marker_xn_vy.outform.header = "UpperLeftVyNodeIndex_xn"
    data.marker_yn_vy.outform.header = "UpperLeftVyNodeIndex_yn"
    data.marker_dx_vy.outform.header = "dx_upperleft_vy"
    data.marker_dy_vy.outform.header = "dy_upperleft_vy"
    return nothing
end

end # module
