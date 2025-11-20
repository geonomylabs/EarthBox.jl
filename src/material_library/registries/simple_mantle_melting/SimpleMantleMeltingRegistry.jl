module SimpleMantleMeltingRegistry

import ..MaterialCollectionManager: MaterialCollection

function get_material_collection()::MaterialCollection
    materials = get_material_names()
    module_dir_path = dirname(@__FILE__)
    path = joinpath(module_dir_path, "simple_mantle_melting.yml")
    return MaterialCollection(materials, path)
end

function get_material_names()::NamedTuple
    return (
        sticky_air = "StickyAir",
        sticky_water = "StickyWater",
        sediment = "Sediment",
        asthenosphere = "Asthenosphere",
        gabbro = "Gabbro",
        gabbroic_magma = "GabbroicMagma"
    )
end

end # module 