module Materials

using EarthBox

const MATERIAL_COLLECTION = MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1

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
        # Sedimentary Basin
        Int16(3) => MaterialDictType(
            "mat_name" => mat_names.clastic_sediment,
            "mat_type" => types.sediment,
            "mat_domain" => domains.sedimentary_basin,
            "red_fraction" => 228.0/255.0,
            "green_fraction" => 148.0/255,
            "blue_fraction" => 73.0/255.0,
        ),
        # Felsic Continental Crust
        Int16(4) => MaterialDictType(
            "mat_name" => mat_names.felsic_continental_crust,
            "mat_type" => types.felsic_continental_crust_fertile,
            "mat_domain" => domains.upper_continental_crust,
            "red_fraction" => 255/255,
            "green_fraction" => 153/255,
            "blue_fraction" => 153/255,
        ),
        Int16(5) => MaterialDictType(
            "mat_name" => mat_names.felsic_continental_crust_strong_zone,
            "mat_type" => types.felsic_continental_crust_fertile,
            "mat_domain" => domains.upper_continental_crust_strong_zone,
            "red_fraction" => 255/255,
            "green_fraction" => 163/255,
            "blue_fraction" => 163/255,
        ),
        # Mafic Continental Crust
        Int16(6) => MaterialDictType(
            "mat_name" => mat_names.mafic_continental_crust,
            "mat_type" => types.mafic_continental_crust_fertile,
            "mat_domain" => domains.lower_continental_crust,
            "red_fraction" => 255/255,
            "green_fraction" => 200/255,
            "blue_fraction" => 200/255,
        ),
        Int16(7) => MaterialDictType(
            "mat_name" => mat_names.mafic_continental_crust_strong_zone,
            "mat_type" => types.mafic_continental_crust_fertile,
            "mat_domain" => domains.lower_continental_crust_strong_zone,
            "red_fraction" => 255/255,
            "green_fraction" => 210/255,
            "blue_fraction" => 210/255,
        ),
        # Continental Mantle Lithosphere
        Int16(8) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_continental_lithosphere,
            "mat_type" => types.ultramafic_mantle_fertile,
            "mat_domain" => domains.mantle_lithosphere,
            "red_fraction" => 0.0,
            "green_fraction" => 153/255,
            "blue_fraction" => 153/255,
        ),
        Int16(9) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_continental_lithosphere_strong_zone,
            "mat_type" => types.ultramafic_mantle_fertile,
            "mat_domain" => domains.lithospheric_mantle_strong_zone,
            "red_fraction" => 0.0,
            "green_fraction" => 163/255,
            "blue_fraction" => 163/255,
        ),
        # Asthenosphere
        Int16(10) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type" => types.ultramafic_mantle_fertile,
            "mat_domain" => domains.asthenosphere,
            "red_fraction" => 0.0,
            "green_fraction" => 200/255,
            "blue_fraction" => 153/255,
        ),
        # Partially Molten Asthenosphere
        Int16(11) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type" => types.ultramafic_mantle_partially_molten,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 255/255,
            "green_fraction" => 150/255,
            "blue_fraction" => 0.0,
        ),
        # Refractory Asthenosphere
        Int16(12) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type" => types.ultramafic_mantle_refactory,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 0.0,
            "green_fraction" => 255/255,
            "blue_fraction" => 0.0,
        ),
        # Solidified Gabbro
        Int16(13) => MaterialDictType(
            "mat_name" => mat_names.oceanic_gabbroic_crust,
            "mat_type" => types.solidified_gabbro,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 227/255,
            "green_fraction" => 227/255,
            "blue_fraction" => 227/255,
        ),
        # Partially Molten Gabbro
        Int16(14) => MaterialDictType(
            "mat_name" => mat_names.oceanic_gabbroic_crust,
            "mat_type" => types.solidified_gabbro_partially_molten,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 255/255,
            "green_fraction" => 125/255,
            "blue_fraction" => 255/255,
        ),
        # Extracted Gabbroic Magma
        Int16(15) => MaterialDictType(
            "mat_name" => mat_names.gabbroic_magma,
            "mat_type" => types.extracted_gabbroic_magma,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 255/255,
            "green_fraction" => 100/255,
            "blue_fraction" => 100/255,
        ),
        # Layered Solidified Gabbroic Crust
        Int16(16) => MaterialDictType(
            "mat_name" => mat_names.layered_gabbroic_crust,
            "mat_type" => types.solidified_layered_gabbro,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 175/255,
            "green_fraction" => 175/255,
            "blue_fraction" => 175/255,
        ),
        # Layered Partially Molten Gabbroic Crust
        Int16(17) => MaterialDictType(
            "mat_name" => mat_names.layered_gabbroic_crust,
            "mat_type" => types.solidified_layered_gabbro_partially_molten,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 255/255,
            "green_fraction" => 0/255,
            "blue_fraction" => 255/255,
        ),
        # Layered Gabbroic Magma
        Int16(18) => MaterialDictType(
            "mat_name" => mat_names.layered_gabbroic_magma,
            "mat_type" => types.extracted_layered_gabbroic_magma,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 255/255,
            "green_fraction" => 0/255,
            "blue_fraction" => 0/255,
        ),
        # Extruded Gabbroic Magma (Lava)
        Int16(19) => MaterialDictType(
            "mat_name" => mat_names.gabbroic_magma,
            "mat_type" => types.extruded_gabbroic_magma,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 255/255,
            "green_fraction" => 255/255,
            "blue_fraction" => 0/255,
        ),
        # Solidified Extruded Gabbroic Magma (Basalt)
        Int16(20) => MaterialDictType(
            "mat_name" => mat_names.oceanic_gabbroic_crust,
            "mat_type" => types.solidified_basalt,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 168/255,
            "green_fraction" => 45/255,
            "blue_fraction" => 0/255,
        ),
    )
end

end # module
