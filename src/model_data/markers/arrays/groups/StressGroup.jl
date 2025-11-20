"""
Module for marker stress arrays.

Provides data structures for storing deviatoric stress components at marker 
locations.
"""
module StressGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "stress"

const ADATA = get_eb_arrays()

"""
    Stress <: AbstractArrayGroup

Array group for marker stress components.

# Fields
- `marker_sxx::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_sxx.description)
- `marker_sxy::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_sxy.description)

# Nested Dot Access
- `marker_sxx = $(ROOT_NAME).$(GRP_NAME).marker_sxx.array`
- `marker_sxy = $(ROOT_NAME).$(GRP_NAME).marker_sxy.array`

# Constructor
    Stress(marknum::Int)

Initializes marker stress arrays with zero values.

## Arguments
- `marknum::Int`: Number of markers

# Returns
- A new Stress struct with the given marker parameters.

"""
mutable struct Stress <: AbstractArrayGroup
    marker_sxx::MarkerArrayFloat1DState{Float64}
    marker_sxy::MarkerArrayFloat1DState{Float64}
end

function Stress(marknum::Int)::Stress
    data = Stress(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),     # array
            ADATA.marker_sxx.name,       # name
            ADATA.marker_sxx.units,      # units
            ADATA.marker_sxx.description # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),     # array
            ADATA.marker_sxy.name,       # name
            ADATA.marker_sxy.units,      # units
            ADATA.marker_sxy.description # description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(stress::Stress)::Nothing
    stress.marker_sxx.outform.fac1 = 1e-6
    stress.marker_sxx.outform.units = "MPa"
    stress.marker_sxx.outform.header = "StressXX_(MPa)"
    stress.marker_sxy.outform.fac1 = 1e-6
    stress.marker_sxy.outform.units = "MPa"
    stress.marker_sxy.outform.header = "StressXY_(MPa)"
    return nothing
end

end # module
