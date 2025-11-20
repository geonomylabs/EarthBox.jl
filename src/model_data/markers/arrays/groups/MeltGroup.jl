"""
Module for marker melt fraction arrays.

Provides data structures for storing melt-related properties of markers 
including total melt fraction, extracted melt, and extractable melt.
"""
module MeltGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "melt"

const ADATA = get_eb_arrays()

"""
    Melt <: AbstractArrayGroup

Array group for marker melt fractions.

# Fields
- `marker_meltfrac::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_meltfrac.description)
- `marker_extracted_meltfrac::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_extracted_meltfrac.description)
- `marker_extractable_meltfrac::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_extractable_meltfrac.description)

# Nested Dot Access
- `marker_meltfrac = $(ROOT_NAME).$(GRP_NAME).marker_meltfrac.array`
- `marker_extracted_meltfrac = $(ROOT_NAME).$(GRP_NAME).marker_extracted_meltfrac.array`
- `marker_extractable_meltfrac = $(ROOT_NAME).$(GRP_NAME).marker_extractable_meltfrac.array`

# Constructor
    Melt(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers

# Returns
- A new Melt struct with the given marker parameters.

"""
mutable struct Melt <: AbstractArrayGroup
    marker_meltfrac::MarkerArrayFloat1DState{Float64}
    marker_extracted_meltfrac::MarkerArrayFloat1DState{Float64}
    marker_extractable_meltfrac::MarkerArrayFloat1DState{Float64}
end

function Melt(marknum::Int)::Melt
    data = Melt(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_meltfrac.name,           # name
            ADATA.marker_meltfrac.units,          # units
            ADATA.marker_meltfrac.description     # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                      # array
            ADATA.marker_extracted_meltfrac.name,        # name
            ADATA.marker_extracted_meltfrac.units,       # units
            ADATA.marker_extracted_meltfrac.description  # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                        # array
            ADATA.marker_extractable_meltfrac.name,        # name
            ADATA.marker_extractable_meltfrac.units,       # units
            ADATA.marker_extractable_meltfrac.description  # description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(data::Melt)::Nothing
    data.marker_meltfrac.outform.header = "MeltFraction"
    data.marker_extracted_meltfrac.outform.header = "ExtractedMeltFrac"
    data.marker_extractable_meltfrac.outform.header = "ExtractableMeltFrac"
    return nothing
end

end # module
