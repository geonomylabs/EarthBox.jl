module SeafloorSpreadingRegistry

import ..MaterialCollectionManager: MaterialCollection

function get_material_collection()::MaterialCollection
    materials = get_material_names()
    module_dir_path = dirname(@__FILE__)
    path = joinpath(module_dir_path, "seafloor_spreading.yml")
    return MaterialCollection(materials, path)
end

function get_material_names()::NamedTuple
    return (
        sfs_sticky_air = "SFSStickyAir",
        sfs_sticky_water = "SFSStickyWater",
        sfs_sediment_wet_quartzite_gt95 = "SFSSedimentWetQuartziteGT95",
        sfs_peridotite_dry_kk08 = "SFSPeridotiteDryKK08",
        sfs_gabbro_wet_anorthite_rd09 = "SFSGabbroWetAnorthiteRD09",
        sfs_gabbroic_magma = "SFSGabbroicMagma",
        sfs_layered_gabbro_wet_anorthite_rd09 = "SFSLayeredGabbroWetAnorthiteRD09",
        sfs_layered_gabbroic_magma = "SFSLayeredGabbroicMagma",
        sfs_serpentinized_peridotite_dry_kk08 = "SFSSerpentinizedPeridotiteDryKK08",
        sfs_upper_continental_crust = "SFSUpperContinentalCrust",
        sfs_lower_continental_crust = "SFSLowerContinentalCrust",
    )
end

end # module 