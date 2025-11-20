module ViscoElastoPlasticRegistry

import ..MaterialCollectionManager: MaterialCollection

function get_material_collection()::MaterialCollection
    materials = get_material_names()
    module_dir_path = dirname(@__FILE__)
    path = joinpath(module_dir_path, "viscoelastoplastic.yml")
    return MaterialCollection(materials, path)
end

function get_material_names()::NamedTuple
    return (
        viscoelastoplastic_sticky_air = "ViscoelastoplasticStickyAir",
        viscoelastoplastic_sticky_water = "ViscoelastoplasticStickyWater",
        viscoelastoplastic_sediment = "ViscoelastoplasticSediment",
        viscoelastoplastic_upper_crust = "ViscoelastoplasticUpperCrust",
        viscoelastoplastic_lower_crust = "ViscoelastoplasticLowerCrust",
        viscoelastoplastic_mantle = "ViscoelastoplasticMantle",
        viscoelastoplastic_weak_fault_crust = "ViscoelastoplasticWeakFaultCrust",
        viscoelastoplastic_weak_fault_mantle = "ViscoelastoplasticWeakFaultMantle"
    )
end

end # module 