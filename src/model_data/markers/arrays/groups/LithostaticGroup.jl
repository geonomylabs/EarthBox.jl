"""
Module for lithostatic-pressure column-filter scratch buffers.

Holds three pre-allocated marker-sized scratch buffers used by
`LithostaticPressure.filter_markers_for_column` to pack column-matching
markers in front. The function then copies the prefix into 3 small tight
return vectors.

Each buffer is single-use within `filter_markers_for_column`, never
shared with any other call site, and fully overwritten each call. Buffers
are passed as **explicit positional arguments** (no
`Union{T, Nothing}=nothing` keyword-arg dispatch hazard).
"""
module LithostaticGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "lithostatic"

const ADATA = get_eb_arrays()

"""
    Lithostatic <: AbstractArrayGroup

Array group holding the three marker-sized scratch buffers used by
`LithostaticPressure.filter_markers_for_column`.

# Fields
- `marker_x_filter_scratch::`[`MarkerArrayFloat1DState`](@ref) Float64
- `marker_y_filter_scratch::`[`MarkerArrayFloat1DState`](@ref) Float64
- `marker_rho_filter_scratch::`[`MarkerArrayFloat1DState`](@ref) Float64

# Constructor
    Lithostatic(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers (each buffer is sized to this length).
"""
mutable struct Lithostatic <: AbstractArrayGroup
    marker_x_filter_scratch::MarkerArrayFloat1DState{Float64}
    marker_y_filter_scratch::MarkerArrayFloat1DState{Float64}
    marker_rho_filter_scratch::MarkerArrayFloat1DState{Float64}
end

function Lithostatic(marknum::Int)::Lithostatic
    return Lithostatic(
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
