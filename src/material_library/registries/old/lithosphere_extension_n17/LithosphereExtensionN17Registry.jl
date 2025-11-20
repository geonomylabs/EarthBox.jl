module LithosphereExtensionN17Registry

"""
    MaterialsRegistry

Registry for materials in library.
"""
@enum MaterialsRegistry begin
    StickyAir = 0
    StickyWater = 1
    Sediment_Wet_Quartzite_GT95 = 2
    Felsic_Crust_Wet_Quartzite_GT95 = 3
    Felsic_Crust_Wet_Quartzite_GT95_Strong_Zone = 4
    Mafic_Crust_Wet_Anorthite_R06 = 5
    Mafic_Crust_Wet_Anorthite_R06_Strong_Zone = 6
    Dry_Olivine_HK03_Lithosphere = 7
    Dry_Olivine_HK03_Lithosphere_Strong_Zone = 8
    Dry_Olivine_HK03_Asthenosphere = 9
    Gabbro_Wet_Anorthite_R06 = 10
    GabbroicMagma = 11
    Felsic_Magma = 12
end

"""
    MaterialNames

Registry for material names in library.

This struct provides easy access to material names in the library.
"""
struct MaterialNames
    sticky_air::String
    sticky_water::String
    sediment_wet_quartzite_gt95::String
    felsic_crust_wet_quartzite_gt95::String
    felsic_crust_wet_quartzite_gt95_strong_zone::String
    mafic_crust_wet_anorthite_r06::String
    mafic_crust_wet_anorthite_r06_strong_zone::String
    dry_olivine_hk03_lithosphere::String
    dry_olivine_hk03_lithosphere_strong_zone::String
    dry_olivine_hk03_asthenosphere::String
    gabbro_wet_anorthite_r06::String
    gabbroic_magma::String
    felsic_magma::String
end

"""
    LithosphereExtensionN17

Material library with registered material names and standard paths.
"""
struct LithosphereExtensionN17
    materials::MaterialNames
    path::String
end

function MaterialNames()
    return MaterialNames(
        string(StickyAir),
        string(StickyWater),
        string(Sediment_Wet_Quartzite_GT95),
        string(Felsic_Crust_Wet_Quartzite_GT95),
        string(Felsic_Crust_Wet_Quartzite_GT95_Strong_Zone),
        string(Mafic_Crust_Wet_Anorthite_R06),
        string(Mafic_Crust_Wet_Anorthite_R06_Strong_Zone),
        string(Dry_Olivine_HK03_Lithosphere),
        string(Dry_Olivine_HK03_Lithosphere_Strong_Zone),
        string(Dry_Olivine_HK03_Asthenosphere),
        string(Gabbro_Wet_Anorthite_R06),
        string(GabbroicMagma),
        string(Felsic_Magma)
    )
end

function LithosphereExtensionN17()
    materials = MaterialNames()
    module_dir_path = dirname(@__FILE__)
    path = joinpath(module_dir_path, "lithosphere_extension_n17.yml")
    return LithosphereExtensionN17(materials, path)
end

end # module 