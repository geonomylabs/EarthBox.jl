module LithosphericDeformationRegistry

import ..MaterialCollectionManager: MaterialCollection

"""
    LithosphericDeformationCollections

The lithospheric deformation collections are a collection of material collections
for the lithospheric deformation models.

See yamal files in src/material_library/registries/lithospheric_deformation.

# Fields
- `lithospheric_deformation_eb1`::[`MaterialCollection`](@ref):
    - Material names and collection file path for the lithospheric deformation model *eb1*
       with moderate crustal flow strength (GT95; R06; no scaling factors), mantle 
       flow laws from HK03, moderate plastic yield strength (cohesion 20 to 10 Pa from plastic strain 0 to 0.1), 
       moderate radiogenic heat production (upper crust 1.8e-6 microW/m^3; lower crust 0.5e-6 microW/m^3).
- `lithospheric_deformation_naliboff17`::[`MaterialCollection`](@ref):
    - Material names and collection file path for the lithospheric deformation model *naliboff17corrected* 
       based on Naliboff et al. (2017) but with temperature-pressure dependent 
       rock properties and melt extraction and corrected pre-exponential term for 
       wet anorthite dislocation creep that includes the water fugacity term.
- `lithospheric_deformation_brune14`::[`MaterialCollection`](@ref):
    - Material names and collection file path for the lithospheric deformation model *brune14corrected*
       based on Brune et al. (2014) but with temperature-pressure dependent 
       rock properties and melt extraction and corrected pre-exponential terms for 
       dislocation and diffusion creep.
"""
struct LithosphericDeformationCollections
    lithospheric_deformation_eb1::MaterialCollection
    lithospheric_deformation_naliboff17::MaterialCollection
    lithospheric_deformation_brune14::MaterialCollection
end

function LithosphericDeformationCollections()
    materials = get_material_names()
    return LithosphericDeformationCollections(
        MaterialCollection(materials, get_collection_path("lithospheric_deformation", "eb1")),
        MaterialCollection(materials, get_collection_path("lithospheric_deformation", "naliboff17")),
        MaterialCollection(materials, get_collection_path("lithospheric_deformation", "brune14")),
    )
end

function get_material_names()::NamedTuple
    return (
        sticky_air = "StickyAir",
        sticky_water = "StickyWater",
        clastic_sediment = "ClasticSediment",
        salt = "Salt",
        felsic_continental_crust = "FelsicContinentalCrust",
        felsic_continental_crust_weak = "FelsicContinentalCrustWeak",
        felsic_continental_crust_strong_zone = "FelsicContinentalCrustStrongZone",
        mafic_continental_crust = "MaficContinentalCrust",
        mafic_continental_crust_strong_zone = "MaficContinentalCrustStrongZone",
        ultramafic_continental_lithosphere = "UltramaficContinentalLithosphere",
        ultramafic_continental_lithosphere_weak = "UltramaficContinentalLithosphereWeak",
        ultramafic_continental_lithosphere_strong_zone = "UltramaficContinentalLithosphereStrongZone",
        ultramafic_asthenosphere_dry_fertile = "UltramaficAsthenosphereDryFertile",
        ultramafic_asthenosphere_wet_fertile = "UltramaficAsthenosphereWetFertile",
        oceanic_gabbroic_crust = "OceanicGabbroicCrust",
        gabbroic_magma = "GabbroicMagma",
        layered_gabbroic_crust = "LayeredGabbroicCrust",
        layered_gabbroic_magma = "LayeredGabbroicMagma",
        serpentinized_peridotite = "SerpentinizedPeridotite",
    )
end

function get_collection_path(group_name::String, model_flag::String)::String
    lib_file_name = string(group_name, "_", model_flag, ".yml")
    module_dir_path = dirname(@__FILE__)
    return joinpath(module_dir_path, lib_file_name)
end

end # module 