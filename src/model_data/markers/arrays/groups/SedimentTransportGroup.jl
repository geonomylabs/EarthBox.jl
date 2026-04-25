"""
Module for sediment-transport marker advection scratch buffers.

Holds pre-allocated marker-sized buffers used by the four advection
functions in `MarkerAdvection.jl`:
- `calculate_sediment_compaction_displacement_factors!`
- `calculate_sediment_marker_displacement!`
- `calculate_sticky_compaction_displacement_factors!`
- `calculate_sticky_marker_displacement!`

The two buffers are reused across sediment and sticky paths within a
single `advect_markers_using_compaction` call. The reuse is safe because:
(a) within a single advect call, the factors buffer is written then
read (no cross-mutation); the displacement buffer likewise; (b) the
sediment path fully completes before the sticky path starts; and (c) each
writer calls `fill!(buffer, 0.0)` first so values never leak between
phases.

Buffers are passed as **explicit positional arguments** to the helpers
(no `Union{T, Nothing}=nothing` keyword), preserving Julia
specialization on the hot path.
"""
module SedimentTransportGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "sediment_transport"

const ADATA = get_eb_arrays()

"""
    SedimentTransport <: AbstractArrayGroup

Array group holding scratch buffers reused by sediment-transport marker
advection routines.

# Fields
- `marker_displacement_factors_buffer::`[`MarkerArrayFloat1DState`](@ref) Float64
- `marker_displacement_buffer::`[`MarkerArrayFloat1DState`](@ref) Float64

# Constructor
    SedimentTransport(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers (each buffer is sized to this length).
"""
mutable struct SedimentTransport <: AbstractArrayGroup
    marker_displacement_factors_buffer::MarkerArrayFloat1DState{Float64}
    marker_displacement_buffer::MarkerArrayFloat1DState{Float64}
end

function SedimentTransport(marknum::Int)::SedimentTransport
    return SedimentTransport(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.sediment_transport_marker_displacement_factors_buffer.name,
            ADATA.sediment_transport_marker_displacement_factors_buffer.units,
            ADATA.sediment_transport_marker_displacement_factors_buffer.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.sediment_transport_marker_displacement_buffer.name,
            ADATA.sediment_transport_marker_displacement_buffer.units,
            ADATA.sediment_transport_marker_displacement_buffer.description
        )
    )
end

end # module
