"""
Module for marker stratigraphy arrays.

Provides data structures for storing stratigraphic information of markers 
including age.
"""
module StratigraphyGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "stratigraphy"

const ADATA = get_eb_arrays()

"""
    Stratigraphy <: AbstractArrayGroup

Array group for marker stratigraphic properties.

# Fields
- `marker_age::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_age.description)

# Nested Dot Access
- `marker_age = $(ROOT_NAME).$(GRP_NAME).marker_age.array`

# Constructor
    Stratigraphy(marknum::Int)

Initializes marker age array with zero values.

## Arguments
- `marknum::Int`: Number of markers

# Returns
- A new Stratigraphy struct with the given marker parameters.

"""
mutable struct Stratigraphy <: AbstractArrayGroup
    marker_age::MarkerArrayFloat1DState{Float64}
end

function Stratigraphy(marknum::Int)::Stratigraphy
    data = Stratigraphy(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),     # array
            ADATA.marker_age.name,       # name
            ADATA.marker_age.units,      # units
            ADATA.marker_age.description # description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(data::Stratigraphy)::Nothing
    data.marker_age.outform.header = "Age_(Myr)"
    return nothing
end

end # module
