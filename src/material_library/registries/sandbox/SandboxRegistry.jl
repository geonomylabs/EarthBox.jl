module SandboxRegistry

import ..MaterialCollectionManager: MaterialCollection

function get_material_collection()::MaterialCollection
    materials = get_sandbox_names()
    module_dir_path = dirname(@__FILE__)
    path = joinpath(module_dir_path, "sandbox.yml")
    return MaterialCollection(materials, path)
end

function get_sandbox_names()::NamedTuple
    return (
        sandbox_sticky_air = "SandboxStickyAir",
        sandbox_sand_layer = "SandboxSandLayer",
        sandbox_microbeads = "SandboxMicrobeads",
        sandbox_mobile_wall = "SandboxMobileWall",
        sandbox_pdms = "SandboxPDMS",
    )
end

end # module 