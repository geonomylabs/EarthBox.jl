"""
Module for serpentinization scratch buffers.

Single marker-sized buffer reused by
`Serpentinization.calculate_marker_serpentinization` to avoid a marknum-
sized `Vector{Float64}` allocation per timestep.
"""
module SerpentinizationGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "serpentinization"

const ADATA = get_eb_arrays()

"""
    Serpentinization <: AbstractArrayGroup

Array group holding scratch buffers reused by serpentinization routines.

# Fields
- `marker_serpentinization_increment_buffer::`[`MarkerArrayFloat1DState`](@ref) Float64

# Constructor
    Serpentinization(marknum::Int)
"""
mutable struct Serpentinization <: AbstractArrayGroup
    marker_serpentinization_increment_buffer::MarkerArrayFloat1DState{Float64}
end

function Serpentinization(marknum::Int)::Serpentinization
    return Serpentinization(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.marker_serpentinization_increment_buffer.name,
            ADATA.marker_serpentinization_increment_buffer.units,
            ADATA.marker_serpentinization_increment_buffer.description
        )
    )
end

end # module
