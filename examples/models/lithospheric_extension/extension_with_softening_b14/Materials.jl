module Materials

using EarthBox

const MATERIAL_COLLECTION = MaterialLibrary().lithospheric_deformation.lithospheric_deformation_brune14

function get_materials_input_dict()::MaterialsDictType
    mat_names = MATERIAL_COLLECTION.materials
    types = MaterialTypesRegistry()
    domains = MaterialDomainsRegistry()

    return MaterialsDictType(
        # Sticky Domain
        Int16(1) => MaterialDictType(
            "mat_name" => mat_names.sticky_air,
            "mat_type" => types.sticky_air,
            "mat_domain" => domains.atmosphere,
            "red_fraction" => 255/255,
            "green_fraction" => 255/255,
            "blue_fraction" => 255/255,
        ),
        Int16(2) => MaterialDictType(
            "mat_name" => mat_names.sticky_water,
            "mat_type" => types.sticky_water,
            "mat_domain" => domains.ocean,
            "red_fraction" => 0/255,
            "green_fraction" => 255/255,
            "blue_fraction" => 255/255,
        ),
        # Felsic Continental Crust
        Int16(3) => MaterialDictType(
            "mat_name" => mat_names.felsic_continental_crust,
            "mat_type" => types.general,
            "mat_domain" => domains.upper_continental_crust,
            "red_fraction" => 255/255,
            "green_fraction" => 153/255,
            "blue_fraction" => 153/255,
        ),
        # Mafic Continental Crust
        Int16(4) => MaterialDictType(
            "mat_name" => mat_names.mafic_continental_crust,
            "mat_type" => types.general,
            "mat_domain" => domains.lower_continental_crust,
            "red_fraction" => 255/255,
            "green_fraction" => 200/255,
            "blue_fraction" => 200/255,
        ),
        # Continental Mantle Lithosphere
        Int16(5) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_continental_lithosphere,
            "mat_type" => types.general,
            "mat_domain" => domains.mantle_lithosphere,
            "red_fraction" => 0.0,
            "green_fraction" => 153/255,
            "blue_fraction" => 153/255,
        ),
        # Asthenosphere
        Int16(6) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type" => types.general,
            "mat_domain" => domains.asthenosphere,
            "red_fraction" => 0.0,
            "green_fraction" => 200/255,
            "blue_fraction" => 153/255,
        )
    )
end

end # module
