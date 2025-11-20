module MarkerWeightsGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.interpolation.arrays"
const GRP_NAME = "marker_weights"

const ADATA = get_eb_arrays()

"""
    MarkerWeights <: AbstractArrayGroup

Array group containing individual marker weight arrays for interpolation.

# Fields
- `marker_wtforULnode::`[`MarkerArrayFloat1DState`](@ref): $(ADATA.marker_wtforULnode.description)
- `marker_wtforLLnode::`[`MarkerArrayFloat1DState`](@ref): $(ADATA.marker_wtforLLnode.description)
- `marker_wtforURnode::`[`MarkerArrayFloat1DState`](@ref): $(ADATA.marker_wtforURnode.description)
- `marker_wtforLRnode::`[`MarkerArrayFloat1DState`](@ref): $(ADATA.marker_wtforLRnode.description)
- `marker_wtforCnode::`[`MarkerArrayFloat1DState`](@ref): $(ADATA.marker_wtforCnode.description)
- `marker_wtforULnodeVy::`[`MarkerArrayFloat1DState`](@ref): $(ADATA.marker_wtforULnodeVy.description)
- `marker_wtforLLnodeVy::`[`MarkerArrayFloat1DState`](@ref): $(ADATA.marker_wtforLLnodeVy.description)
- `marker_wtforURnodeVy::`[`MarkerArrayFloat1DState`](@ref): $(ADATA.marker_wtforURnodeVy.description)
- `marker_wtforLRnodeVy::`[`MarkerArrayFloat1DState`](@ref): $(ADATA.marker_wtforLRnodeVy.description)

# Nested Dot Access
- `marker_wtforULnode = $(ROOT_NAME).$(GRP_NAME).marker_wtforULnode.array`
- `marker_wtforLLnode = $(ROOT_NAME).$(GRP_NAME).marker_wtforLLnode.array`
- `marker_wtforURnode = $(ROOT_NAME).$(GRP_NAME).marker_wtforURnode.array`
- `marker_wtforLRnode = $(ROOT_NAME).$(GRP_NAME).marker_wtforLRnode.array`
- `marker_wtforCnode = $(ROOT_NAME).$(GRP_NAME).marker_wtforCnode.array`
- `marker_wtforULnodeVy = $(ROOT_NAME).$(GRP_NAME).marker_wtforULnodeVy.array`
- `marker_wtforLLnodeVy = $(ROOT_NAME).$(GRP_NAME).marker_wtforLLnodeVy.array`
- `marker_wtforURnodeVy = $(ROOT_NAME).$(GRP_NAME).marker_wtforURnodeVy.array`
- `marker_wtforLRnodeVy = $(ROOT_NAME).$(GRP_NAME).marker_wtforLRnodeVy.array`

# Constructor
    MarkerWeights(marknum::Int)

# Arguments
- `marknum::Int`: Number of markers in the model

# Returns
- `MarkerWeights`: New MarkerWeights array group with initialized arrays

"""
mutable struct MarkerWeights <: AbstractArrayGroup
    marker_wtforULnode::MarkerArrayFloat1DState{Float64}
    marker_wtforLLnode::MarkerArrayFloat1DState{Float64}
    marker_wtforURnode::MarkerArrayFloat1DState{Float64}
    marker_wtforLRnode::MarkerArrayFloat1DState{Float64}
    marker_wtforCnode::MarkerArrayFloat1DState{Float64}
    marker_wtforULnodeVy::MarkerArrayFloat1DState{Float64}
    marker_wtforLLnodeVy::MarkerArrayFloat1DState{Float64}
    marker_wtforURnodeVy::MarkerArrayFloat1DState{Float64}
    marker_wtforLRnodeVy::MarkerArrayFloat1DState{Float64}
end

function MarkerWeights(marknum::Int)::MarkerWeights
    return MarkerWeights(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                    # array
            ADATA.marker_wtforULnode.name,              # name
            ADATA.marker_wtforULnode.units,             # units
            ADATA.marker_wtforULnode.description        # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                    # array
            ADATA.marker_wtforLLnode.name,              # name
            ADATA.marker_wtforLLnode.units,             # units
            ADATA.marker_wtforLLnode.description        # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                    # array
            ADATA.marker_wtforURnode.name,              # name
            ADATA.marker_wtforURnode.units,             # units
            ADATA.marker_wtforURnode.description        # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                    # array
            ADATA.marker_wtforLRnode.name,              # name
            ADATA.marker_wtforLRnode.units,             # units
            ADATA.marker_wtforLRnode.description        # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                    # array
            ADATA.marker_wtforCnode.name,               # name
            ADATA.marker_wtforCnode.units,              # units
            ADATA.marker_wtforCnode.description         # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                    # array
            ADATA.marker_wtforULnodeVy.name,            # name
            ADATA.marker_wtforULnodeVy.units,           # units
            ADATA.marker_wtforULnodeVy.description      # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                    # array
            ADATA.marker_wtforLLnodeVy.name,            # name
            ADATA.marker_wtforLLnodeVy.units,           # units
            ADATA.marker_wtforLLnodeVy.description      # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                    # array
            ADATA.marker_wtforURnodeVy.name,            # name
            ADATA.marker_wtforURnodeVy.units,           # units
            ADATA.marker_wtforURnodeVy.description      # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                    # array
            ADATA.marker_wtforLRnodeVy.name,            # name
            ADATA.marker_wtforLRnodeVy.units,           # units
            ADATA.marker_wtforLRnodeVy.description      # description
        )
    )
end

end # module 