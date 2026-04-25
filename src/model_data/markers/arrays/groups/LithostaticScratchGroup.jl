"""
Module for lithostatic-pressure column-filter scratch buffers.

Holds three pre-allocated marker-sized buffers used by
`LithostaticPressure.filter_markers_for_column` to avoid 3 marknum-sized
`Vector{Float64}` allocations per sealevel update.
"""
module LithostaticScratchGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "lithostatic_scratch"

const ADATA = get_eb_arrays()

"""
    LithostaticScratch <: AbstractArrayGroup

Array group of marker-sized scratch buffers reused by the column-filter
inside `LithostaticPressure.filter_markers_for_column`.

# Fields
- `marker_x_filter_scratch::`[`MarkerArrayFloat1DState`](@ref) Float64
- `marker_y_filter_scratch::`[`MarkerArrayFloat1DState`](@ref) Float64
- `marker_rho_filter_scratch::`[`MarkerArrayFloat1DState`](@ref) Float64

# Constructor
    LithostaticScratch(marknum::Int)
"""
mutable struct LithostaticScratch <: AbstractArrayGroup
    marker_x_filter_scratch::MarkerArrayFloat1DState{Float64}
    marker_y_filter_scratch::MarkerArrayFloat1DState{Float64}
    marker_rho_filter_scratch::MarkerArrayFloat1DState{Float64}
end

function LithostaticScratch(marknum::Int)::LithostaticScratch
    return LithostaticScratch(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.marker_x_filter_scratch.name,
            ADATA.marker_x_filter_scratch.units,
            ADATA.marker_x_filter_scratch.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.marker_y_filter_scratch.name,
            ADATA.marker_y_filter_scratch.units,
            ADATA.marker_y_filter_scratch.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.marker_rho_filter_scratch.name,
            ADATA.marker_rho_filter_scratch.units,
            ADATA.marker_rho_filter_scratch.description
        )
    )
end

end # module
