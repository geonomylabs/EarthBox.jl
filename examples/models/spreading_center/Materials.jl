module Materials

using EarthBox

const MATERIAL_COLLECTION = MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1

function get_materials_input_dict()::MaterialsDictType
    mat_names = MATERIAL_COLLECTION.materials
    types = MaterialTypesRegistry()
    domains = MaterialDomainsRegistry()

    return MaterialsDictType(
        # Sticky Air
        Int16(1) => MaterialDictType(
            "mat_name" => mat_names.sticky_air,
            "mat_type" => types.sticky_air,
            "mat_domain" => domains.atmosphere,
            "red_fraction" => 215/255,
            "green_fraction" => 215/255,
            "blue_fraction" => 215/255,
        ),
        # Sticky Water
        Int16(2) => MaterialDictType(
            "mat_name" => mat_names.sticky_water,
            "mat_type" => types.sticky_water,
            "mat_domain" => domains.ocean,
            "red_fraction" => 0/255,
            "green_fraction" => 255/255,
            "blue_fraction" => 255/255,
        ),
        # Clastic Sediment
        Int16(3) => MaterialDictType(
            "mat_name" => mat_names.clastic_sediment,
            "mat_type" => types.sediment,
            "mat_domain" => domains.sedimentary_basin,
            "red_fraction" => 204/255,
            "green_fraction" => 102/255,
            "blue_fraction" => 0/255,
        ),
        # Asthenosphere
        Int16(4) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type" => types.ultramafic_mantle_fertile,
            "mat_domain" => domains.asthenosphere,
            "red_fraction" => 0/255,
            "green_fraction" => 204/255,
            "blue_fraction" => 102/255,
        ),
        # Partially Molten Asthenosphere
        Int16(5) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type" => types.ultramafic_mantle_partially_molten,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 255/255,
            "green_fraction" => 153/255,
            "blue_fraction" => 51/255,
        ),
        # Refractory Asthenosphere
        Int16(6) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type" => types.ultramafic_mantle_refactory,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 0/255,
            "green_fraction" => 255/255,
            "blue_fraction" => 0/255,
        ),
        # Solidified Gabbro
        Int16(7) => MaterialDictType(
            "mat_name" => mat_names.oceanic_gabbroic_crust,
            "mat_type" => types.solidified_gabbro,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 128/255,
            "green_fraction" => 128/255,
            "blue_fraction" => 128/255,
        ),
        # Partially Molten Gabbro
        Int16(8) => MaterialDictType(
            "mat_name" => mat_names.oceanic_gabbroic_crust,
            "mat_type" => types.solidified_gabbro_partially_molten,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 127/255,
            "green_fraction" => 0/255,
            "blue_fraction" => 255/255,
        ),
        # Extracted Gabbroic Magma
        Int16(9) => MaterialDictType(
            "mat_name" => mat_names.gabbroic_magma,
            "mat_type" => types.extracted_gabbroic_magma,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 189/255,
            "green_fraction" => 0/255,
            "blue_fraction" => 0/255,
        ),
        # Layered Solidified Gabbroic Crust
        Int16(10) => MaterialDictType(
            "mat_name" => mat_names.layered_gabbroic_crust,
            "mat_type" => types.solidified_layered_gabbro,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 175/255,
            "green_fraction" => 175/255,
            "blue_fraction" => 175/255,
        ),
        # Layered Partially Molten Gabbroic Crust
        Int16(11) => MaterialDictType(
            "mat_name" => mat_names.layered_gabbroic_crust,
            "mat_type" => types.solidified_layered_gabbro_partially_molten,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 255/255,
            "green_fraction" => 0/255,
            "blue_fraction" => 255/255,
        ),
        # Layered Gabbroic Magma
        Int16(12) => MaterialDictType(
            "mat_name" => mat_names.layered_gabbroic_magma,
            "mat_type" => types.extracted_layered_gabbroic_magma,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 255/255,
            "green_fraction" => 0/255,
            "blue_fraction" => 0/255,
        ),
        # Extruded Gabbroic Magma (Lava)
        Int16(13) => MaterialDictType(
            "mat_name" => mat_names.gabbroic_magma,
            "mat_type" => types.extruded_gabbroic_magma,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 255/255,
            "green_fraction" => 255/255,
            "blue_fraction" => 0/255,
        ),
        # Solidified Extruded Gabbroic Magma (Basalt)
        Int16(14) => MaterialDictType(
            "mat_name" => mat_names.oceanic_gabbroic_crust,
            "mat_type" => types.solidified_basalt,
            "mat_domain" => domains.general_domain,
            "red_fraction" => 0/255,
            "green_fraction" => 0/255,
            "blue_fraction" => 0/255,
        ),
    )
end

end # module
