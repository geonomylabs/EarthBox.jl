module Gerya2019Registry

import ..MaterialCollectionManager: MaterialCollection

function get_material_collection()::MaterialCollection
    materials = get_material_names()
    module_dir_path = dirname(@__FILE__)
    path = joinpath(module_dir_path, "gerya2019.yml")
    return MaterialCollection(materials, path)
end

function get_material_names()::NamedTuple
    return (
        sticky_air = "StickyAir",
        sticky_water = "StickyWater",
        sediment = "Sediment",
        basalt = "Basalt",
        gabbro = "Gabbro",
        mantle = "Mantle",
        hydrated_mantle = "HydratedMantle",
    )
end

end # module 