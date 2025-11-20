module BenchmarksRegistry

import EarthBox.TupleTools: merge_named_tuples
import ..MaterialCollectionManager: MaterialCollection

function get_material_collection()::MaterialCollection
    materials = get_material_names()
    module_dir_path = dirname(@__FILE__)
    path = joinpath(module_dir_path, "benchmarks.yml")
    return MaterialCollection(materials, path)
end

function get_material_names()::NamedTuple
    return merge_named_tuples(
        get_viscoelastic_names(),
        get_solid_body_rotation_names(),
        get_rayleigh_taylor_names(),
        get_isoviscous_fluid_names(),
        get_viscoelastic_slab_names(),
        get_kaus10_names(),
        get_nonnewtonian_names(),
        get_constant_viscosity_channel_names(),
        get_temp_dependent_viscosity_couette_flow_names(),
        get_flextest_names(),
        get_seafloor_spreading_names(),
        get_extension_names(),
        get_forced_adiabatic_gradient_names(),
        get_viscoelastic_extension_names(),
        get_simple_sedimentation_names(),
    )
end

function get_viscoelastic_names()::NamedTuple
    return (
        viscoelastic_material = "Viscoelastic_Material",
    )
end

function get_solid_body_rotation_names()::NamedTuple
    return (
        material_for_solid_body_rotation = "MaterialForSolidBodyRotation",
    )
end

function get_rayleigh_taylor_names()::NamedTuple
    return (
        layer1_rayleigh_taylor = "Layer1RayleighTaylor",
        layer2_rayleigh_taylor = "Layer2RayleighTaylor",
    )
end

function get_isoviscous_fluid_names()::NamedTuple
    return (
        isoviscous_fluid_blankenbach89 = "IsoviscousFluidBlankenbach89",
        isoviscous_fluid_1_around_elastic_slab = "IsoviscousFluid1AroundElasticSlab",
        isoviscous_fluid_2_around_elastic_slab = "IsoviscousFluid2AroundElasticSlab",
    )
end

function get_viscoelastic_slab_names()::NamedTuple
    return (
        viscoelastic_slab_material_1 = "ViscoelasticSlabMaterial1",
        viscoelastic_slab_material_2 = "ViscoelasticSlabMaterial2",
    )
end

function get_kaus10_names()::NamedTuple
    return (
        sticky_air_kaus10 = "StickyAirKaus10",
        crust_kaus10 = "CrustKaus10",
        weak_seed_kaus10 = "WeakSeedKaus10",
    )
end

function get_nonnewtonian_names()::NamedTuple
    return (
        nonnewtonian_material_channel_flow_benchmark = "NonNewtonianMaterialChannelFlowBenchmark",
    )
end

function get_constant_viscosity_channel_names()::NamedTuple
    return (
        constant_viscosity_channel = "ConstantViscosityChannel",
        constant_viscosity_channel_variable_thermal_conductivity = "ConstantViscosityChannelVariableThermalConductivity",
    )
end

function get_temp_dependent_viscosity_couette_flow_names()::NamedTuple

    return (
        temp_dependent_viscosity_couette_flow_benchmark = "TempDependentViscosityCouetteFlowBenchmark",
    )
end

function get_flextest_names()::NamedTuple
    return (
        flexteststickywater = "FlexTestStickyWater",
        flextestisoviscouscrust = "FlexTestIsoviscousCrust",
        flextestisoviscousmantle = "FlexTestIsoviscousMantle",
    )
end

function get_seafloor_spreading_names()::NamedTuple
    return (
        sfs_sticky_air = "SFSStickyAir",
        sfs_sticky_water = "SFSStickyWater",
        sfs_clastic_sediment = "SFSClasticSediment",
        sfs_peridotite_dry_fertile = "SFSPeridotiteDryFertile",
        sfs_gabbroic_crust = "SFSGabbroicCrust",
        sfs_gabbroic_magma = "SFSGabbroicMagma",
        sfs_layered_gabbroic_crust = "SFSLayeredGabbroicCrust",
        sfs_layered_gabbroic_magma = "SFSLayeredGabbroicMagma",
        sfs_upper_continental_crust = "SFSUpperContinentalCrust",
        sfs_lower_continental_crust = "SFSLowerContinentalCrust",
    )
end

function get_extension_names()::NamedTuple
    return (
        extension_sticky_air = "ExtensionStickyAir",
        extension_sticky_water = "ExtensionStickyWater",
        extension_sediment_wet_quartzite_gt95 = "ExtensionSedimentWetQuartziteGT95",
        extension_sediment_wet_quartzite_r95 = "ExtensionSedimentWetQuartziteR95",
        extension_felsic_crust_wet_quartzite_gt95 = "ExtensionFelsicCrustWetQuartziteGT95",
        extension_felsic_crust_wet_quartzite_gt95_strong_zone = "ExtensionFelsicCrustWetQuartziteGT95StrongZone",
        extension_mafic_crust_wet_anorthite_rd09 = "ExtensionMaficCrustWetAnorthiteRD09",
        extension_mafic_crust_wet_anorthite_rd09_strong_zone = "ExtensionMaficCrustWetAnorthiteRD09StrongZone",
        extension_dry_olivine_kk08_strong_zone = "ExtensionDryOlivineKK08StrongZone",
        extension_peridotite_dry_kk08 = "ExtensionPeridotiteDryKK08",
        extension_dry_olivine_kk08_weak_seed = "ExtensionDryOlivineKK08WeakSeed",
        extension_wet_olivine_kk08 = "ExtensionWetOlivineKK08",
        extension_felsic_magma = "ExtensionFelsicMagma",
        extension_gabbroic_magma = "ExtensionGabbroicMagma",
        extension_gabbro_wet_anorthite_rd09 = "ExtensionGabbroWetAnorthiteRD09",
    )
end

function get_forced_adiabatic_gradient_names()::NamedTuple
    return (
        forced_adiabatic_gradient_sticky_air = "ForcedAdiabaticGradientStickyAir",
        forced_adiabatic_gradient_dry_olivine_hk03_lithosphere = "ForcedAdiabaticGradientDryOlivineHK03Lithosphere",
        forced_adiabatic_gradient_dry_olivine_hk03_lithosphere_strong_zone = "ForcedAdiabaticGradientDryOlivineHK03LithosphereStrongZone",
        forced_adiabatic_gradient_dry_olivine_hk03_asthenosphere = "ForcedAdiabaticGradientDryOlivineHK03Asthenosphere",
    )
end

function get_viscoelastic_extension_names()::NamedTuple
    return (
        viscoelastic_extension_sticky_air = "ViscoelasticExtensionStickyAir",
        viscoelastic_extension_sticky_water = "ViscoelasticExtensionStickyWater",
        viscoelastic_extension_sediment = "ViscoelasticExtensionSediment",
        viscoelastic_extension_upper_crust = "ViscoelasticExtensionUpperCrust",
        viscoelastic_extension_lower_crust = "ViscoelasticExtensionLowerCrust",
        viscoelastic_extension_upper_crust_strong = "ViscoelasticExtensionUpperCrustStrong",
        viscoelastic_extension_lower_crust_strong = "ViscoelasticExtensionLowerCrustStrong",
        viscoelastic_extension_mantle = "ViscoelasticExtensionMantle",
    )
end

function get_simple_sedimentation_names()::NamedTuple
    return (
        simple_sedimentation_sticky_air = "SimpleSedimentationStickyAir",
        simple_sedimentation_sticky_water = "SimpleSedimentationStickyWater",
        simple_sedimentation_sediment_wet_quartzite_gt95 = "SimpleSedimentationSedimentWetQuartzite_GT95",
        simple_sedimentation_peridotite_dry_kk08 = "SimpleSedimentationPeridotiteDryKK08",
        simple_sedimentation_salt = "SimpleSedimentationSalt",
    )
end

end # module 