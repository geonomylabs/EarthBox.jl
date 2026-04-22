"""
Module for marker advection velocity and spin arrays.

Provides storage for Runge-Kutta-interpolated velocity and spin values at
marker locations. These arrays are filled each timestep by the advection
velocity interpolator and then consumed by the marker-location and
stress-rotation updates.
"""
module AdvectionGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "advection"

const ADATA = get_eb_arrays()

"""
    Advection <: AbstractArrayGroup

Array group for per-marker Runge-Kutta-interpolated velocity and spin. The
arrays are pre-allocated and re-written every advection step; consumers must
gate reads by the same `inside_flags` used by the producer.

# Fields
- `marker_vx::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_vx.description)
- `marker_vy::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_vy.description)
- `marker_spin::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_spin.description)

# Nested Dot Access
- `marker_vx = $(ROOT_NAME).$(GRP_NAME).marker_vx.array`
- `marker_vy = $(ROOT_NAME).$(GRP_NAME).marker_vy.array`
- `marker_spin = $(ROOT_NAME).$(GRP_NAME).marker_spin.array`

# Constructor
    Advection(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers.
"""
mutable struct Advection <: AbstractArrayGroup
    marker_vx::MarkerArrayFloat1DState{Float64}
    marker_vy::MarkerArrayFloat1DState{Float64}
    marker_spin::MarkerArrayFloat1DState{Float64}
end

function Advection(marknum::Int)::Advection
    data = Advection(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.marker_vx.name,
            ADATA.marker_vx.units,
            ADATA.marker_vx.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.marker_vy.name,
            ADATA.marker_vy.units,
            ADATA.marker_vy.description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),
            ADATA.marker_spin.name,
            ADATA.marker_spin.units,
            ADATA.marker_spin.description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(data::Advection)::Nothing
    data.marker_vx.outform.header = "MarkerVx"
    data.marker_vy.outform.header = "MarkerVy"
    data.marker_spin.outform.header = "MarkerSpin"
    return nothing
end

end # module
