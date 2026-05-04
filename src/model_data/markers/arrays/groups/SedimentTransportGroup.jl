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

Both buffers are constructed at length 0 and lazily resized to `marknum`
on the first call to `advect_markers_using_compaction` via
[`ensure_sediment_transport_buffers!`](@ref). They are flagged with
`ibackup=false` so the JLD2 backup loader does not try to restore them
from saved state — they are scratch and re-derived on each advect call.
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
    SedimentTransport()

The marker-length buffers start empty and are sized on first use by
[`ensure_sediment_transport_buffers!`](@ref).
"""
mutable struct SedimentTransport <: AbstractArrayGroup
    marker_displacement_factors_buffer::MarkerArrayFloat1DState{Float64}
    marker_displacement_buffer::MarkerArrayFloat1DState{Float64}
end

function SedimentTransport()::SedimentTransport
    return SedimentTransport(
        MarkerArrayFloat1DState(
            Float64[],
            ADATA.sediment_transport_marker_displacement_factors_buffer.name,
            ADATA.sediment_transport_marker_displacement_factors_buffer.units,
            ADATA.sediment_transport_marker_displacement_factors_buffer.description;
            ibackup=false
        ),
        MarkerArrayFloat1DState(
            Float64[],
            ADATA.sediment_transport_marker_displacement_buffer.name,
            ADATA.sediment_transport_marker_displacement_buffer.units,
            ADATA.sediment_transport_marker_displacement_buffer.description;
            ibackup=false
        )
    )
end

"""
    ensure_sediment_transport_buffers!(st::SedimentTransport, marknum::Int)

Idempotently size the marker-length advection scratch buffers to
`marknum`. Safe to call every advect — does nothing when sizes already
match. Called at the top of
`MarkerAdvection.advect_markers_using_compaction`.
"""
function ensure_sediment_transport_buffers!(
    st::SedimentTransport, marknum::Int
)::Nothing
    _ensure_buffer!(st.marker_displacement_factors_buffer.array, marknum)
    _ensure_buffer!(st.marker_displacement_buffer.array, marknum)
    return nothing
end

function _ensure_buffer!(buf::Vector{Float64}, target::Int)::Nothing
    current = length(buf)
    if current == target
        return nothing
    end
    resize!(buf, target)
    if target > current
        @inbounds for i in (current + 1):target
            buf[i] = 0.0
        end
    end
    return nothing
end

end # module
