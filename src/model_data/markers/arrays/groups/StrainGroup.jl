"""
Module for marker strain and strain rate arrays.

Provides data structures for storing strain-related properties of markers 
including accumulated strain, strain rates, and melt damage.
"""
module StrainGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "strain"

const ADATA = get_eb_arrays()

"""
    Strain <: AbstractArrayGroup

Array group for marker strain and strain rate.

# Fields
- `marker_GII::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_GII.description)
- `marker_strain_plastic::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_strain_plastic.description)
- `marker_exx::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_exx.description)
- `marker_exy::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_exy.description)
- `marker_sr_ratio::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_sr_ratio.description)
- `marker_strain_rate_plastic::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_strain_rate_plastic.description)
- `marker_melt_damage::`[`MarkerArrayFloat1DState`](@ref) Float32: $(ADATA.marker_melt_damage.description)

# Nested Dot Access
- `marker_GII = $(ROOT_NAME).$(GRP_NAME).marker_GII.array`
- `marker_strain_plastic = $(ROOT_NAME).$(GRP_NAME).marker_strain_plastic.array`
- `marker_exx = $(ROOT_NAME).$(GRP_NAME).marker_exx.array`
- `marker_exy = $(ROOT_NAME).$(GRP_NAME).marker_exy.array`
- `marker_sr_ratio = $(ROOT_NAME).$(GRP_NAME).marker_sr_ratio.array`
- `marker_strain_rate_plastic = $(ROOT_NAME).$(GRP_NAME).marker_strain_rate_plastic.array`
- `marker_melt_damage = $(ROOT_NAME).$(GRP_NAME).marker_melt_damage.array`

# Constructor
    Strain(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers

# Returns
- A new Strain struct with the given marker parameters.

"""
mutable struct Strain <: AbstractArrayGroup
    marker_GII::MarkerArrayFloat1DState{Float64}
    marker_strain_plastic::MarkerArrayFloat1DState{Float64}
    marker_exx::MarkerArrayFloat1DState{Float64}
    marker_exy::MarkerArrayFloat1DState{Float64}
    marker_sr_ratio::MarkerArrayFloat1DState{Float64}
    marker_strain_rate_plastic::MarkerArrayFloat1DState{Float64}
    marker_melt_damage::MarkerArrayFloat1DState{Float64}
end

function Strain(marknum::Int)::Strain
    data = Strain(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_GII.name,                # name
            ADATA.marker_GII.units,               # units
            ADATA.marker_GII.description          # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                      # array
            ADATA.marker_strain_plastic.name,             # name
            ADATA.marker_strain_plastic.units,            # units
            ADATA.marker_strain_plastic.description       # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_exx.name,                # name
            ADATA.marker_exx.units,               # units
            ADATA.marker_exx.description          # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_exy.name,                # name
            ADATA.marker_exy.units,               # units
            ADATA.marker_exy.description          # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_sr_ratio.name,           # name
            ADATA.marker_sr_ratio.units,          # units
            ADATA.marker_sr_ratio.description     # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                      # array
            ADATA.marker_strain_rate_plastic.name,        # name
            ADATA.marker_strain_rate_plastic.units,       # units
            ADATA.marker_strain_rate_plastic.description  # description
        ),
        MarkerArrayFloat1DState(
            ones(Float64, marknum),                       # array
            ADATA.marker_melt_damage.name,                # name
            ADATA.marker_melt_damage.units,               # units
            ADATA.marker_melt_damage.description          # description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(data::Strain)::Nothing
    data.marker_GII.outform.header = "Strain"
    data.marker_strain_plastic.outform.header = "PlasticStrain"
    data.marker_sr_ratio.outform.header = "EIImarker/EIIgrid"
    data.marker_exx.outform.header = "EXX_(1/s)"
    data.marker_exx.outform.log10 = true
    data.marker_exy.outform.header = "EXY_(1/s)"
    data.marker_exy.outform.log10 = true
    data.marker_strain_rate_plastic.outform.header = "PlasticStrainRate(1/s)"
    data.marker_strain_rate_plastic.outform.log10 = true
    data.marker_melt_damage.outform.header = "MeltDamage"
    return nothing
end

end # module
