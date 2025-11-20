"""
Module for marker thermal property arrays.

Provides data structures for storing thermal properties of markers including
temperature, heat capacity, thermal conductivity, and adiabatic heating.
"""
module ThermalGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "thermal"

const ADATA = get_eb_arrays()

"""
    Thermal <: AbstractArrayGroup

Array group for marker thermal properties.

# Fields
- `marker_TK::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_TK.description)
- `marker_rhocp::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_rhocp.description)
- `marker_kt::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_kt.description)
- `marker_ha::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_ha.description)

# Nested Dot Access
- `marker_TK = $(ROOT_NAME).$(GRP_NAME).marker_TK.array`
- `marker_rhocp = $(ROOT_NAME).$(GRP_NAME).marker_rhocp.array`
- `marker_kt = $(ROOT_NAME).$(GRP_NAME).marker_kt.array`
- `marker_ha = $(ROOT_NAME).$(GRP_NAME).marker_ha.array`

# Constructor
    Thermal(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers

# Returns
- A new Thermal struct with the given marker parameters.

"""
mutable struct Thermal <: AbstractArrayGroup
    marker_TK::MarkerArrayFloat1DState{Float64}
    marker_rhocp::MarkerArrayFloat1DState{Float64}
    marker_kt::MarkerArrayFloat1DState{Float64}
    marker_ha::MarkerArrayFloat1DState{Float64}
end

function Thermal(marknum::Int)::Thermal
    data = Thermal(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),     # array
            ADATA.marker_TK.name,        # name
            ADATA.marker_TK.units,       # units
            ADATA.marker_TK.description  # description
        ),
        MarkerArrayFloat1DState(
            fill(3300.0 * 1100.0, marknum), # array
            ADATA.marker_rhocp.name,         # name
            ADATA.marker_rhocp.units,        # units
            ADATA.marker_rhocp.description   # description
        ),
        MarkerArrayFloat1DState(
            fill(3.0, marknum),          # array
            ADATA.marker_kt.name,        # name
            ADATA.marker_kt.units,       # units
            ADATA.marker_kt.description  # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),     # array
            ADATA.marker_ha.name,        # name
            ADATA.marker_ha.units,       # units
            ADATA.marker_ha.description  # description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(thermal::Thermal)::Nothing
    thermal.marker_TK.outform.fac2 = -273.0
    thermal.marker_TK.outform.units = "C"
    thermal.marker_TK.outform.header = "Temperature_(C)"
    thermal.marker_rhocp.outform.header = "RhoCp_(J/K/m^3)"
    thermal.marker_kt.outform.header = "ThermalCond_(W/m/K)"
    thermal.marker_ha.outform.header = "AdiabaticHeatingTerm_(alphaT)"
    return nothing
end

end # module
