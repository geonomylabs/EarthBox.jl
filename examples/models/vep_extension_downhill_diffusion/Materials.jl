module Materials

using EarthBox

const MATERIAL_COLLECTION = MaterialLibrary().viscoelastoplastic

function get_materials_input_dict()::MaterialsDictType
    mat_names = MATERIAL_COLLECTION.materials
    types = MaterialTypesRegistry()
    domains = MaterialDomainsRegistry()

    return MaterialsDictType(
        # Sticky Air
        Int16(1) => MaterialDictType(
            "mat_name" => mat_names.viscoelastoplastic_sticky_air,
            "mat_type" => types.sticky_air,
            "mat_domain" => domains.atmosphere,
            "red_fraction" => 1.0,
            "green_fraction" => 1.0,
            "blue_fraction" => 1.0,
        ),
        # Sticky Water
        Int16(2) => MaterialDictType(
            "mat_name" => mat_names.viscoelastoplastic_sticky_water,
            "mat_type" => types.sticky_water,
            "mat_domain" => "Ocean",
            "red_fraction" => 0.0,
            "green_fraction" => 1.0,
            "blue_fraction" => 1.0,
        ),
        # Sediment
        Int16(3) => MaterialDictType(
            "mat_name" => mat_names.viscoelastoplastic_sediment,
            "mat_type" => types.sediment,
            "mat_domain" => domains.sedimentary_basin,
            "red_fraction" => 228.0/255.0,
            "green_fraction" => 148.0/255.0,
            "blue_fraction" => 73.0/255.0,
        ),
        # Upper Crust
        Int16(4) => MaterialDictType(
            "mat_name" => mat_names.viscoelastoplastic_upper_crust,
            "mat_type" => types.general,
            "mat_domain" => domains.upper_continental_crust,
            "red_fraction" => 1.0,
            "green_fraction" => 0.6,
            "blue_fraction" => 0.6,
        ),
        # Lower Crust
        Int16(5) => MaterialDictType(
            "mat_name" => mat_names.viscoelastoplastic_lower_crust,
            "mat_type" => types.general,
            "mat_domain" => domains.lower_continental_crust,
            "red_fraction" => 1.0,
            "green_fraction" => 0.8627451,
            "blue_fraction" => 0.8627451,
        ),
        # Upper Mantle Lithosphere
        Int16(6) => MaterialDictType(
            "mat_name" => mat_names.viscoelastoplastic_mantle,
            "mat_type" => types.general,
            "mat_domain" => domains.mantle_lithosphere,
            "red_fraction" => 0.0,
            "green_fraction" => 153/255,
            "blue_fraction" => 153/255,
        ),
        # Asthenosphere
        Int16(7) => MaterialDictType(
            "mat_name" => mat_names.viscoelastoplastic_mantle,
            "mat_type" => types.general,
            "mat_domain" => domains.asthenosphere,
            "red_fraction" => 0.0,
            "green_fraction" => 200/255,
            "blue_fraction" => 153/255,
        ),
        # Weak Fault Crust
        Int16(8) => MaterialDictType(
            "mat_name" => mat_names.viscoelastoplastic_weak_fault_crust,
            "mat_type" => types.felsic_continental_crust_fertile,
            "mat_domain" => domains.weak_crustal_fault_zone,
            "red_fraction" => 1.0,
            "green_fraction" => 128.0/255.0,
            "blue_fraction" => 0.0,
        ),
        # Weak Fault Mantle
        Int16(9) => MaterialDictType(
            "mat_name" => mat_names.viscoelastoplastic_weak_fault_mantle,
            "mat_type" => types.general,
            "mat_domain" => domains.weak_mantle_fault_zone,
            "red_fraction" => 177.0/255.0,
            "green_fraction" => 1.0,
            "blue_fraction" => 0.0,
        )
    )
end

end # module
