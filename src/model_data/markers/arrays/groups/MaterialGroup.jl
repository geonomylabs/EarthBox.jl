"""
Module for marker material property arrays.

Provides data structures for storing material properties of markers including
material ID, density, porosity, and serpentinization.
"""
module MaterialGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "material"

const ADATA = get_eb_arrays()

"""
    Material <: AbstractArrayGroup

Array group for marker material properties.

# Fields
- `marker_matid::`[`MarkerArrayInt1DState`](@ref) Int16: $(ADATA.marker_matid.description)
- `marker_rho::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_rho.description)
- `marker_porosity_initial::`[`MarkerArrayFloat1DState`](@ref) Float32: $(ADATA.marker_porosity_initial.description)
- `marker_decay_depth::`[`MarkerArrayFloat1DState`](@ref) Float32: $(ADATA.marker_decay_depth.description)
- `marker_max_burial_depth::`[`MarkerArrayFloat1DState`](@ref) Float32: $(ADATA.marker_max_burial_depth.description)
- `marker_serpentinization::`[`MarkerArrayFloat1DState`](@ref) Float32: $(ADATA.marker_serpentinization.description)
- `marker_serpentinization_heat_production::`[`MarkerArrayFloat1DState`](@ref) Float32: $(ADATA.marker_serpentinization_heat_production.description)

# Nested Dot Access
- `marker_matid = $(ROOT_NAME).$(GRP_NAME).marker_matid.array`
- `marker_rho = $(ROOT_NAME).$(GRP_NAME).marker_rho.array`
- `marker_porosity_initial = $(ROOT_NAME).$(GRP_NAME).marker_porosity_initial.array`
- `marker_decay_depth = $(ROOT_NAME).$(GRP_NAME).marker_decay_depth.array`
- `marker_max_burial_depth = $(ROOT_NAME).$(GRP_NAME).marker_max_burial_depth.array`
- `marker_serpentinization = $(ROOT_NAME).$(GRP_NAME).marker_serpentinization.array`
- `marker_serpentinization_heat_production = $(ROOT_NAME).$(GRP_NAME).marker_serpentinization_heat_production.array`

# Constructor
    Material(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers

# Returns
- A new Material struct with the given marker parameters.

"""
mutable struct Material <: AbstractArrayGroup
    marker_matid::MarkerArrayInt1DState{Int16}
    marker_rho::MarkerArrayFloat1DState{Float64}
    marker_porosity_initial::MarkerArrayFloat1DState{Float32}
    marker_decay_depth::MarkerArrayFloat1DState{Float32}
    marker_max_burial_depth::MarkerArrayFloat1DState{Float32}
    marker_serpentinization::MarkerArrayFloat1DState{Float64}
    marker_serpentinization_heat_production::MarkerArrayFloat1DState{Float64}
end

function Material(marknum::Int)::Material
    data = Material(
        MarkerArrayInt1DState(
            zeros(Int16, marknum),                    # array
            ADATA.marker_matid.name,                  # name
            ADATA.marker_matid.units,                 # units
            ADATA.marker_matid.description            # description
        ),
        MarkerArrayFloat1DState(
            fill(3300.0, marknum),                    # array
            ADATA.marker_rho.name,                    # name
            ADATA.marker_rho.units,                   # units
            ADATA.marker_rho.description              # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float32, marknum),                  # array
            ADATA.marker_porosity_initial.name,       # name
            ADATA.marker_porosity_initial.units,      # units
            ADATA.marker_porosity_initial.description # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float32, marknum),                  # array
            ADATA.marker_decay_depth.name,            # name
            ADATA.marker_decay_depth.units,           # units
            ADATA.marker_decay_depth.description      # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float32, marknum),                  # array
            ADATA.marker_max_burial_depth.name,       # name
            ADATA.marker_max_burial_depth.units,      # units
            ADATA.marker_max_burial_depth.description # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                  # array
            ADATA.marker_serpentinization.name,       # name
            ADATA.marker_serpentinization.units,      # units
            ADATA.marker_serpentinization.description # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),                              # array
            ADATA.marker_serpentinization_heat_production.name,   # name
            ADATA.marker_serpentinization_heat_production.units,  # units
            ADATA.marker_serpentinization_heat_production.description # description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(
    data::Material
)::Nothing
    data.marker_matid.outform.header = "MaterialID"
    data.marker_rho.outform.header = "Density_(kg/m^3)"
    data.marker_porosity_initial.outform.header = "Initial_Porosity"
    data.marker_decay_depth.outform.header = "Decay_Depth_(m)"
    data.marker_max_burial_depth.outform.header = "Max_Burial_Depth_(m)"
    data.marker_serpentinization.outform.header = "Serpentinization_Fraction"
    data.marker_serpentinization_heat_production.outform.header = "Serpentinization_Heat_Production_Heat_Production_(W/m^3)"
    return nothing
end

end # module