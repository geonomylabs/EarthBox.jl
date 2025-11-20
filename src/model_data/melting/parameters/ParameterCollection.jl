module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import .OptionsGroup: Options
import .MeltingMaterialIDsGroup: MeltingMaterialIDs
import .RheologyGroup: Rheology
import .ExtrusionGroup: Extrusion
import .ExtractionGroup: Extraction

"""
    Parameters <: AbstractParameterCollection

Parameter collection for melting and melt extraction processes.

# Fields
- `options::`[`Options`](@ref): Melting model on/off switches
- `melting_material_ids::`[`MeltingMaterialIDs`](@ref): Material IDs for melting processes
- `rheology::`[`Rheology`](@ref): Rheological properties of molten rock
- `extrusion::`[`Extrusion`](@ref): Parameters for melt extrusion
- `extraction::`[`Extraction`](@ref): Parameters for melt extraction

# Constructor
    Parameters()

Create a new Parameters collection with default melting model values.

"""
mutable struct Parameters <: AbstractParameterCollection
    options::Options
    melting_material_ids::MeltingMaterialIDs
    rheology::Rheology
    extrusion::Extrusion
    extraction::Extraction
end

function Parameters()::Parameters
    return Parameters(
        Options(),
        MeltingMaterialIDs(),
        Rheology(),
        Extrusion(),
        Extraction()
    )
end

end # module 