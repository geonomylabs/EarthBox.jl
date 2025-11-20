"""
Module for marker pressure arrays.

Provides data structures for storing pressure values at marker locations.
"""
module PressureGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "pressure"

const ADATA = get_eb_arrays()

"""
    Pressure <: AbstractArrayGroup

Array group for marker pressure values.

# Fields
- `marker_pr::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_pr.description)

# Nested Dot Access
- `marker_pr = $(ROOT_NAME).$(GRP_NAME).marker_pr.array`

# Constructor
    Pressure(marknum::Int)

Initializes marker pressure array with zero values.

## Arguments
- `marknum::Int`: Number of markers
"""
mutable struct Pressure <: AbstractArrayGroup
    marker_pr::MarkerArrayFloat1DState{Float64}
end

function Pressure(marknum::Int)::Pressure
    data = Pressure(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),     # array
            ADATA.marker_pr.name,        # name
            ADATA.marker_pr.units,       # units
            ADATA.marker_pr.description  # description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(data::Pressure)::Nothing
    data.marker_pr.outform.fac1 = 1e-9
    data.marker_pr.outform.units = "GPa"
    data.marker_pr.outform.header = "Pressure_(GPa)"
    return nothing
end

end # module
