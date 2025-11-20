"""
Module for marker rheological property arrays.

Provides data structures for storing rheological properties of markers 
including viscosity, friction, cohesion, and plastic failure indicators.
"""
module RheologyGroup

import EarthBox.Arrays.ArrayRegistry: get_eb_arrays
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.EarthBoxDtypes: AbstractArrayGroup

const ROOT_NAME = "model.markers.arrays"
const GRP_NAME = "rheology"

const ADATA = get_eb_arrays()

"""
    Rheology <: AbstractArrayGroup

Array group for marker rheological properties.

# Fields
- `marker_eta::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_eta.description)
- `marker_fric_ini::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_fric_ini.description)
- `marker_fric::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_fric.description)
- `marker_cohesion::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_cohesion.description)
- `marker_preexp::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_preexp.description)
- `marker_pfailure::`[`MarkerArrayFloat1DState`](@ref) Float32: $(ADATA.marker_pfailure.description)
- `marker_eta_flow::`[`MarkerArrayFloat1DState`](@ref) Float64: $(ADATA.marker_eta_flow.description)
- `marker_dilatation_angle::`[`MarkerArrayFloat1DState`](@ref) Float32: $(ADATA.marker_dilatation_angle.description)

# Nested Dot Access
- `marker_eta = $(ROOT_NAME).$(GRP_NAME).marker_eta.array`
- `marker_fric_ini = $(ROOT_NAME).$(GRP_NAME).marker_fric_ini.array`
- `marker_fric = $(ROOT_NAME).$(GRP_NAME).marker_fric.array`
- `marker_cohesion = $(ROOT_NAME).$(GRP_NAME).marker_cohesion.array`
- `marker_preexp = $(ROOT_NAME).$(GRP_NAME).marker_preexp.array`
- `marker_pfailure = $(ROOT_NAME).$(GRP_NAME).marker_pfailure.array`
- `marker_eta_flow = $(ROOT_NAME).$(GRP_NAME).marker_eta_flow.array`
- `marker_dilatation_angle = $(ROOT_NAME).$(GRP_NAME).marker_dilatation_angle.array`

# Constructor
    Rheology(marknum::Int)

## Arguments
- `marknum::Int`: Number of markers

# Returns
- A new Rheology struct with the given marker parameters.

"""
mutable struct Rheology <: AbstractArrayGroup
    marker_eta::MarkerArrayFloat1DState{Float64}
    marker_fric_ini::MarkerArrayFloat1DState{Float64}
    marker_fric::MarkerArrayFloat1DState{Float64}
    marker_cohesion::MarkerArrayFloat1DState{Float64}
    marker_preexp::MarkerArrayFloat1DState{Float64}
    marker_pfailure::MarkerArrayFloat1DState{Float32}
    marker_eta_flow::MarkerArrayFloat1DState{Float64}
    marker_dilatation_angle::MarkerArrayFloat1DState{Float32}
end

function Rheology(marknum::Int)::Rheology
    data = Rheology(
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_eta.name,                # name
            ADATA.marker_eta.units,               # units
            ADATA.marker_eta.description          # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_fric_ini.name,           # name
            ADATA.marker_fric_ini.units,          # units
            ADATA.marker_fric_ini.description     # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_fric.name,               # name
            ADATA.marker_fric.units,              # units
            ADATA.marker_fric.description         # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_cohesion.name,           # name
            ADATA.marker_cohesion.units,          # units
            ADATA.marker_cohesion.description     # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_preexp.name,             # name
            ADATA.marker_preexp.units,            # units
            ADATA.marker_preexp.description       # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float32, marknum),              # array
            ADATA.marker_pfailure.name,           # name
            ADATA.marker_pfailure.units,          # units
            ADATA.marker_pfailure.description     # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float64, marknum),              # array
            ADATA.marker_eta_flow.name,           # name
            ADATA.marker_eta_flow.units,          # units
            ADATA.marker_eta_flow.description     # description
        ),
        MarkerArrayFloat1DState(
            zeros(Float32, marknum),                      # array
            ADATA.marker_dilatation_angle.name,           # name
            ADATA.marker_dilatation_angle.units,          # units
            ADATA.marker_dilatation_angle.description     # description
        )
    )
    update_output_format!(data)
    return data
end

function update_output_format!(data::Rheology)::Nothing
    data.marker_eta.outform.header = "Viscosity_(Pa.s)"
    data.marker_fric_ini.outform.header = "Initial_Friction_Coeff"
    data.marker_fric.outform.header = "Friction_Coeff"
    data.marker_cohesion.outform.header = "Cohesion_(Pa)"
    data.marker_preexp.outform.header = "Pre-expo_(1/s/MPa^n)"
    data.marker_pfailure.outform.header = "Plastic_Failure_Flag"
    data.marker_eta_flow.outform.header = "Viscosity_FlowLaw_(Pa.s)"
    data.marker_dilatation_angle.outform.header = "Dilatation_Angle_Degrees"
    return nothing
end

end # module
