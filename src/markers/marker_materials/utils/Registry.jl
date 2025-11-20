module Registry

"""
    MaterialTypesRegistry

A registry of material type names. Material types are used to define different 
evolving genetic characteristics of a material including phase transformation 
and other evolving genetic characteristics. This structure provides easy access
to the material type names.

# Fields
- `general::String` = "General".
    - General material type. This is the default material type for all materials.
- `sticky_air::String` = "StickyAir".
    - Sticky air material type used to approximate a stress-free air-rock interface that may involve
      the transformation of rock material to sticky air to approximate the effects of erosion.
- `sticky_water::String` = "StickyWater".
    - Sticky water material type used to approximate a stress-free water-rock interface that may involve
      the transformation of rock material to sticky water to approximate the effects of erosion or
      the transformation of sticky-air to sticky-water to maintain the water-air interface.
- `sediment::String` = "Sediment".
    - Sediment material type formed by erosion and deposition often implemented in the code by 
      transforming sticky air/water to sediment to approximate deposition.
- `felsic_continental_crust_fertile::String` = "FertileFelsicContinentalCrust".
   - Fertile felsic continental crust capable of undergoing partial melting to produce magma.
- `felsic_continental_crust_partially_molten::String` = "PartiallyMoltenFelsicContinentalCrust".
   - Partially molten felsic continental crust containing a mixture of solid and liquid phases.
- `felsic_continental_crust_refactory::String` = "RefactoryFelsicContinentalCrust".
   - Refactory felsic continental crust that has undergone melting and exists under conditions not 
      conducive to further melting.
- `extracted_felsic_magma::String` = "ExtractedFelsicMagma".
    - Melt extracted from the felsic continental crust to form intra-crustal magma bodies.
- `solidified_granite::String` = "SolidifiedGranite".
    - Solidified granite material type formed by the cooling of magma.
- `solidified_rhyolite::String` = "SolidifiedRhyolite".
    - Solidified rhyolite material type formed by the cooling of granitic magma extruded onto the 
      Earth's surface.
- `mafic_continental_crust_fertile::String` = "FertileMaficContinentalCrust".
   - Fertile mafic continental crust capable of undergoing partial melting to produce magma.
- `mafic_continental_crust_partially_molten::String` = "PartiallyMoltenMaficContinentalCrust".
   - Partially molten mafic continental crust containing a mixture of solid and liquid phases.
- `mafic_continental_crust_refactory::String` = "RefactoryMaficContinentalCrust".
   - Refactory mafic continental crust that has undergone melting and exists under conditions not 
      conducive to further melting.
- `ultramafic_mantle_fertile::String` = "FertileUltramaficMantle".
   - Fertile ultramafic mantle capable of undergoing partial melting to produce magma.
- `ultramafic_mantle_partially_molten::String` = "PartiallyMoltenUltramaficMantle".
   - Partially molten ultramafic mantle containing a mixture of solid and liquid phases.
- `ultramafic_mantle_refactory::String` = "RefactoryUltramaficMantle".
   - Refactory ultramafic mantle that has undergone melting and exists under conditions not 
      conducive to further melting.
- `ultramafic_continental_mantle_refractory::String` = "UltramaficContinentalMantleRefractory".
   - Refractory ultramafic mantle that has undergone melting and exists under conditions not 
      conducive to further melting.
- `solidified_gabbro::String` = "SolidifiedGabbro".
    - Solidified gabbro material type formed by the cooling of gabbroic magma or partially molten 
       gabbro.
- `solidified_gabbro_partially_molten::String` = "SolidifiedGabbroPartiallyMolten".
    - Partially molten gabbro material that has undergone partial solidification from gabbroic magma
      and still contains a mixture of solid and liquid phases.
- `extracted_gabbroic_magma::String` = "ExtractedGabbroicMagma".
    - Gabbroic melt extracted from the mantle and instantaneously migrated to 
       to form magma bodies in accumulation zones either at local maxima in the 
       partially molten region or at the base of the crust. 
- `extruded_gabbroic_magma::String` = "ExtrudedGabbroicMagma".
    - Gabbroic melt extracted from the mantle and instantly extruded onto the 
      Earth's surface.
- `solidified_layered_gabbro::String` = "SolidifiedLayeredGabbro".
    - Solidified layered gabbro material type formed by the cooling of layered gabbroic magma or 
       partially molten layered gabbro. Layered gabbro refers to gabbroic material that has undergone
       fractional crystallization that modifies the composition and melting properties.
- `solidified_layered_gabbro_partially_molten::String` = "SolidifiedLayeredGabbroPartiallyMolten".
    - Partially molten layered gabbro material that has undergone partial solidification from layered 
       gabbroic magma and still contains a mixture of solid and liquid phases.
- `extracted_layered_gabbroic_magma::String` = "ExtractedLayeredGabbroicMagma".
    - Layered gabbroic melt extracted from the mantle and instantaneously migrated to 
       to form magma bodies in accumulation zones at the base of the Moho above local
       maxima in the partially molten mantle regions. 
- `solidified_basalt::String` = "SolidifiedBasalt".
    - Solidified basalt material type formed by the cooling of basaltic magma or partially molten 
       basalt.
- `serpentinite::String` = "Serpentinite".
    - Serpentinite material type formed by the alteration of peridotite in the presence of water
      and implemented in the code by transforming ultramafic material to serpentinite.
- `salt::String` = "Salt".
    - Salt material type formed by the evaporation and precipitation and implemented in the code by 
      transforming sticky water to salt.

"""
Base.@kwdef struct MaterialTypesRegistry
    general::String = "General"
    sticky_air::String = "StickyAir"
    sticky_water::String = "StickyWater"
    sediment::String = "Sediment"
    # Felsic Continental Crustal Rocks
    felsic_continental_crust_fertile::String = "FelsicContinentalCrustFertile"
    felsic_continental_crust_partially_molten::String = "FelsicContinentalCrustPartiallyMolten"
    felsic_continental_crust_refactory::String = "FelsicContinentalCrustRefactory"
    extracted_felsic_magma::String = "ExtractedFelsicMagma"
    solidified_granite::String = "SolidifiedGranite"
    solidified_rhyolite::String = "SolidifiedRhyolite"
    # Mafic Continental Crustal Rocks
    mafic_continental_crust_fertile::String = "MaficContinentalCrustFertile"
    mafic_continental_crust_partially_molten::String = "MaficContinentalCrustPartiallyMolten"
    mafic_continental_crust_refactory::String = "MaficContinentalCrustRefactory"
    # Ultramafic Mantle Rocks
    ultramafic_mantle_fertile::String = "UltramaficMantleFertile"
    ultramafic_mantle_partially_molten::String = "UltramaficMantlePartiallyMolten"
    ultramafic_mantle_refactory::String = "UltramaficMantleRefactory"
    ultramafic_continental_mantle_refractory::String = "UltramaficContinentalMantleRefractory"
    # Gabbro
    solidified_gabbro::String = "SolidifiedGabbro"
    solidified_gabbro_partially_molten::String = "SolidifiedGabbroPartiallyMolten"
    extracted_gabbroic_magma::String = "ExtractedGabbroicMagma"
    extruded_gabbroic_magma::String = "ExtrudedGabbroicMagma"
    # LayeredGabbro (Fractionated gabbro in lower oceanic crust)
    solidified_layered_gabbro::String = "SolidifiedLayeredGabbro"
    solidified_layered_gabbro_partially_molten::String = "SolidifiedLayeredGabbroPartiallyMolten"
    extracted_layered_gabbroic_magma::String = "ExtractedLayeredGabbroicMagma"
    # Basalt
    solidified_basalt::String = "SolidifiedBasalt"
    # Serpentinite
    serpentinite::String = "Serpentinite"
    # Salt
    salt::String = "Salt"
end

function get_types_registry()::Vector{String}
    mtr = MaterialTypesRegistry()
    return [
        mtr.general,
        mtr.sticky_air,
        mtr.sticky_water,
        mtr.sediment,
        mtr.felsic_continental_crust_fertile,
        mtr.felsic_continental_crust_partially_molten,
        mtr.felsic_continental_crust_refactory,
        mtr.extracted_felsic_magma,
        mtr.solidified_granite,
        mtr.solidified_rhyolite,
        mtr.mafic_continental_crust_fertile,
        mtr.mafic_continental_crust_partially_molten,
        mtr.mafic_continental_crust_refactory,
        mtr.ultramafic_mantle_fertile,
        mtr.ultramafic_mantle_partially_molten,
        mtr.ultramafic_mantle_refactory,
        mtr.ultramafic_continental_mantle_refractory,
        mtr.solidified_gabbro,
        mtr.solidified_gabbro_partially_molten,
        mtr.extracted_gabbroic_magma,
        mtr.extruded_gabbroic_magma,
        mtr.solidified_layered_gabbro,
        mtr.solidified_layered_gabbro_partially_molten,
        mtr.extracted_layered_gabbroic_magma,
        mtr.solidified_basalt,
        mtr.serpentinite,
        mtr.salt
    ]
end

"""
    MaterialDomainsRegistry

A registry of material domains. Material domains are used to define the initial 
physical domain of a material. Material domains are used by initialization 
functions to assign materials to the correct physical domain.

# Fields
- `general_domain::String` = "GeneralDomain".
    - General material domain. This is the default material domain for all materials.
- `atmosphere::String` = "Atmosphere".
    - Atmosphere domain composed of sticky-air material.
- `ocean::String` = "Ocean".
    - Ocean domain composed of sticky-water material.
- `sedimentary_basin::String` = "SedimentaryBasin".
    - Sedimentary basin domain composed of sediment material.
- `continental_crust::String` = "ContinentalCrust".
    - Continental crust domain composed of felsic rocks.
- `upper_continental_crust::String` = "UpperContinentalCrust".
    - Upper continental crust domain composed of felsic rocks.
- `lower_continental_crust::String` = "LowerContinentalCrust".
    - Lower continental crust domain composed of felsic and/or mafic crustal rocks.
- `exhumed_continental_crust::String` = "ExhumedContinentalCrust".
    - Exhumed continental crust domain composed of felsic and/or mafic rocks that have been
      exhumed from the lower continental crust to the upper continental crust.
- `oceanic_crust::String` = "OceanicCrust".
    - Oceanic crust material domain composed of mafic rocks formed at spreading centers.
- `basaltic_oceanic_crust::String` = "BasalticOceanicCrust".
    - Basaltic oceanic crust domain composed of basaltic rocks formed at spreading centers.
- `gabbroic_oceanic_crust::String` = "GabbroicOceanicCrust".
    - Gabbroic oceanic crust domain composed of gabbroic rocks formed at spreading centers.
- `exhumed_mantle_lithosphere::String` = "ExhumedMantleLithosphere".
    - Exhumed ultramafic mantle lithosphere.
- `mantle_lithosphere::String` = "MantleLithosphere".
    - Ultramafic mantle lithosphere. If this domain is defined in the material model then the domains 
      "UpperMantleLithosphere", "MiddleMantleLithosphere" and "LowerMantleLithosphere"
      will have properties assigned to them automatically from the "MantleLithosphere" domain.
- `weak_mantle_lithosphere::String` = "WeakMantleLithosphere".
    - Weak mantle lithosphere domain with weaker frictional-plastic properties and/or 
      lower relative effective viscosity.
- `upper_mantle_lithosphere::String` = "UpperMantleLithosphere".
    - Upper mantle lithosphere domain composed of ultramafic rocks.
- `middle_mantle_lithosphere::String` = "MiddleMantleLithosphere".
    - Middle mantle lithosphere domain composed of ultramafic rocks.
- `lower_mantle_lithosphere::String` = "LowerMantleLithosphere".
    - Lower mantle lithosphere domain composed of ultramafic rocks.
- `asthenosphere::String` = "Asthenosphere".
    - Asthenosphere domain below the mantle lithosphere composed of ultramafic rocks.
- `mantle::String` = "Mantle".
    - Mantle domain composed of ultramafic rocks.
- `lower_mantle::String` = "LowerMantle".
    - Lower mantle domain composed of ultramafic rocks.
- `mantle_plume::String` = "MantlePlume".
    - Mantle plume domain composed of ultramafic rocks.
- `upper_continental_crust_strong_zone::String` = "UpperContinentalCrustStrongZone".
    - Upper continental crust strong zone with relatively less damage and stronger 
      frictional-plastic properties.
- `lower_continental_crust_strong_zone::String` = "LowerContinentalCrustStrongZone".
    - Lower continental crust strong zone with relatively less damage and stronger 
      frictional-plastic properties.
- `lithospheric_mantle_strong_zone::String` = "LithosphericMantleStrongZone".
    - Lithospheric mantle strong zone with relatively less damage and stronger 
      frictional-plastic properties.
- `weak_seed_crust::String` = "WeakSeedCrust".
    - Weak seed crust domain with weaker frictional-plastic properties and/or 
      lower relative effective viscosity.
- `weak_seed_mantle::String` = "WeakSeedMantle".
    - Weak seed mantle domain with weaker frictional-plastic properties and/or 
      lower relative effective viscosity.
- `weak_crustal_fault_zone::String` = "WeakCrustalFaultZone".
    - Weak crustal fault zone domain with weaker frictional-plastic properties and/or 
      lower relative effective viscosity.
- `weak_mantle_fault_zone::String` = "WeakMantleFaultZone".
    - Weak mantle fault zone domain with weaker frictional-plastic properties and/or 
      lower relative effective viscosity.
- `fluid_around_elastic_slab_a::String` = "FluidAroundElasticSlabA".
    - Fluid domain around elastic slab with unique color for creating a checkerboard
      pattern.
- `fluid_around_elastic_slab_b::String` = "FluidAroundElasticSlabB".
    - Fluid domain around elastic slab with unique color for creating a checkerboard
      pattern.
- `elastic_slab_a::String` = "ElasticSlabA".
    - Elastic slab with unique color for creating a checkerboard pattern.
- `elastic_slab_b::String` = "ElasticSlabB".
    - Elastic slab with unique color for creating a checkerboard pattern.
- `layer1_rayleigh_taylor::String` = "Layer1RayleighTaylor".
    - Layer 1 of the Rayleigh Taylor instability benchmark.
- `layer2_rayleigh_taylor::String` = "Layer2RayleighTaylor".
    - Layer 2 of the Rayleigh Taylor instability benchmark.
- `sand_a::String` = "SandA".
    - Sand with unique color for creating a layered stratigraphic pattern.
- `sand_b::String` = "SandB".
    - Sand domain with unique color for creating a layered stratigraphic pattern.
- `microbeads::String` = "Microbeads".
    - Microbead domain used in sandbox experiments.
- `mobile_wall::String` = "MobileWall".
    - Strong mobile wall domain used in sandbox experiments.
- `pdms_layer::String` = "PDMSLayer".
    - Polydimethylsiloxane (PDMS) layer domain used in sandbox experiments.
- `plate_extension::String` = "PlateExtension".
    - Strong plate domain that extends from mobile wall used in sandbox experiments.
"""
Base.@kwdef struct MaterialDomainsRegistry
    general_domain::String = "GeneralDomain"
    atmosphere::String = "Atmosphere"
    ocean::String = "Ocean"
    sedimentary_basin::String = "SedimentaryBasin"
    continental_crust::String = "ContinentalCrust"
    upper_continental_crust::String = "UpperContinentalCrust"
    lower_continental_crust::String = "LowerContinentalCrust"
    exhumed_continental_crust::String = "ExhumedContinentalCrust"
    oceanic_crust::String = "OceanicCrust"
    basaltic_oceanic_crust::String = "BasalticOceanicCrust"
    gabbroic_oceanic_crust::String = "GabbroicOceanicCrust"
    exhumed_mantle_lithosphere::String = "ExhumedMantleLithosphere"
    mantle_lithosphere::String = "MantleLithosphere"
    weak_mantle_lithosphere::String = "WeakMantleLithosphere"
    upper_mantle_lithosphere::String = "UpperMantleLithosphere"
    middle_mantle_lithosphere::String = "MiddleMantleLithosphere"
    lower_mantle_lithosphere::String = "LowerMantleLithosphere"
    asthenosphere::String = "Asthenosphere"
    mantle::String = "Mantle"
    lower_mantle::String = "LowerMantle"
    mantle_plume::String = "MantlePlume"
    upper_continental_crust_strong_zone::String = "UpperContinentalCrustStrongZone"
    lower_continental_crust_strong_zone::String = "LowerContinentalCrustStrongZone"
    lithospheric_mantle_strong_zone::String = "LithosphericMantleStrongZone"
    weak_seed_crust::String = "WeakSeedCrust"
    weak_seed_mantle::String = "WeakSeedMantle"
    weak_crustal_fault_zone::String = "WeakCrustalFaultZone"
    weak_mantle_fault_zone::String = "WeakMantleFaultZone"
    fluid_around_elastic_slab_a::String = "FluidAroundElasticSlabA"
    fluid_around_elastic_slab_b::String = "FluidAroundElasticSlabB"
    elastic_slab_a::String = "ElasticSlabA"
    elastic_slab_b::String = "ElasticSlabB"
    layer1_rayleigh_taylor::String = "Layer1RayleighTaylor"
    layer2_rayleigh_taylor::String = "Layer2RayleighTaylor"
    sand_a::String = "SandA"
    sand_b::String = "SandB"
    microbeads::String = "Microbeads"
    mobile_wall::String = "MobileWall"
    pdms_layer::String = "PDMSLayer"
    plate_extension::String = "PlateExtension"
end


function get_domains_registry()::Vector{String}
    mdr = MaterialDomainsRegistry()
    return [
        mdr.general_domain,
        mdr.atmosphere,
        mdr.ocean,
        mdr.sedimentary_basin,
        mdr.continental_crust,
        mdr.upper_continental_crust,
        mdr.lower_continental_crust,
        mdr.exhumed_continental_crust,
        mdr.oceanic_crust,
        mdr.basaltic_oceanic_crust,
        mdr.gabbroic_oceanic_crust,
        mdr.exhumed_mantle_lithosphere,
        mdr.mantle_lithosphere,
        mdr.weak_mantle_lithosphere,
        mdr.upper_mantle_lithosphere,
        mdr.middle_mantle_lithosphere,
        mdr.lower_mantle_lithosphere,
        mdr.asthenosphere,
        mdr.mantle,
        mdr.lower_mantle,
        mdr.mantle_plume,
        mdr.upper_continental_crust_strong_zone,
        mdr.lower_continental_crust_strong_zone,
        mdr.lithospheric_mantle_strong_zone,
        mdr.weak_seed_crust,
        mdr.weak_seed_mantle,
        mdr.weak_crustal_fault_zone,
        mdr.weak_mantle_fault_zone,
        mdr.fluid_around_elastic_slab_a,
        mdr.fluid_around_elastic_slab_b,
        mdr.elastic_slab_a,
        mdr.elastic_slab_b,
        mdr.layer1_rayleigh_taylor,
        mdr.layer2_rayleigh_taylor,
        mdr.sand_a,
        mdr.sand_b,
        mdr.microbeads,
        mdr.mobile_wall,
        mdr.pdms_layer,
        mdr.plate_extension
    ]
end

end # module Registry 