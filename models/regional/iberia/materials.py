""" Materials input dict for lithospheric extension and breakup model.
"""
from earthbox import MaterialsLibrary
from earthbox.markers.marker_materials.materials.utils.material_types \
    import MaterialsDictType
from earthbox.markers.marker_materials import MaterialTypesRegistry
from earthbox.markers.marker_materials import MaterialDomainsRegistry


def get_materials_input_dict() -> tuple[MaterialsDictType, str]:
    """ Get materials input dictionary for model and path to materials library.

    This function returns a user defined materials model and path to the
    selected materials library used to build to build this model.

    Returns
    -------
    materials_input_dict : MaterialsDictType
        Materials input dictionary.

    path : str
        Path to material library.

    """

    # Select a materials library
    lib = MaterialsLibrary().lithospheric_deformation.lithospheric_deformation_a1

    mat_names = lib.materials
    types = MaterialTypesRegistry()
    domains = MaterialDomainsRegistry()

    # Define materials input dictionary
    materials_input_dict = {
        # Sticky Domain
        0: {
            "mat_name": mat_names.sticky_air,
            "mat_type": types.sticky_air,
            "mat_domain": domains.atmosphere,
            "red_fraction": 255/255,
            "green_fraction": 255/255,
            "blue_fraction": 255/255
        },
        1: {
            "mat_name": mat_names.sticky_water,
            "mat_type": types.sticky_water,
            "mat_domain": domains.ocean,
            "red_fraction": 0/255,
            "green_fraction": 255/255,
            "blue_fraction": 255/255
        },
        # Sedimentary Basin
        2: {
            "mat_name": mat_names.clastic_sediment,
            "mat_type": types.sediment,
            "mat_domain": domains.sedimentary_basin,
            "red_fraction": 0.89411765,
            "green_fraction": 0.58039216,
            "blue_fraction": 0.28627451
        },
        # Felsic Continental Crust
        3: {
            "mat_name": mat_names.felsic_continental_crust,
            "mat_type": types.felsic_continental_crust_fertile,
            "mat_domain": domains.upper_continental_crust,
            "red_fraction": 255/255,
            "green_fraction": 153/255,
            "blue_fraction": 153/255
        },
        4: {
            "mat_name": mat_names.felsic_continental_crust_strong_zone,
            "mat_type": types.felsic_continental_crust_fertile,
            "mat_domain": domains.upper_continental_crust_strong_zone,
            "red_fraction": 255/255,
            "green_fraction": 163/255,
            "blue_fraction": 163/255
        },
        # Mafic Continental Crust
        5: {
            "mat_name": mat_names.mafic_continental_crust,
            "mat_type": types.mafic_continental_crust_fertile,
            "mat_domain": domains.lower_continental_crust,
            "red_fraction": 255/255,
            "green_fraction": 200/255,
            "blue_fraction": 200/255
        },
        6: {
            "mat_name": mat_names.mafic_continental_crust_strong_zone,
            "mat_type": types.mafic_continental_crust_fertile,
            "mat_domain": domains.lower_continental_crust_strong_zone,
            "red_fraction": 255/255,
            "green_fraction": 210/255,
            "blue_fraction": 210/255
        },
        # Continental Mantle Lithosphere
        7: {
            "mat_name": mat_names.ultramafic_continental_lithosphere,
            "mat_type": types.ultramafic_mantle_fertile,
            "mat_domain": domains.upper_mantle_lithosphere,
            "red_fraction": 0.0,
            "green_fraction": 153/255,
            "blue_fraction": 153/255
        },
        8: {
            "mat_name": mat_names.ultramafic_continental_lithosphere,
            "mat_type": types.ultramafic_mantle_fertile,
            "mat_domain": domains.middle_mantle_lithosphere,
            "red_fraction": 0.0,
            "green_fraction": 153/255,
            "blue_fraction": 153/255
        },
        9: {
            "mat_name": mat_names.ultramafic_continental_lithosphere,
            "mat_type": types.ultramafic_mantle_fertile,
            "mat_domain": domains.lower_mantle_lithosphere,
            "red_fraction": 0.0,
            "green_fraction": 153/255,
            "blue_fraction": 153/255
        },
        10: {
            "mat_name": mat_names.ultramafic_continental_lithosphere_strong_zone,
            "mat_type": types.ultramafic_mantle_fertile,
            "mat_domain": domains.lithospheric_mantle_strong_zone,
            "red_fraction": 0.0,
            "green_fraction": 163/255,
            "blue_fraction": 163/255
        },
        # Asthenosphere
        11: {
            "mat_name": mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type": types.ultramafic_mantle_fertile,
            "mat_domain": domains.asthenosphere,
            "red_fraction": 0.0,
            "green_fraction": 200/255,
            "blue_fraction": 153/255
        },
        # Partially Molten Asthenosphere
        12: {
            "mat_name": mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type": types.ultramafic_mantle_partially_molten,
            "mat_domain": domains.general_domain,
            "red_fraction": 255/255,
            "green_fraction": 150/255,
            "blue_fraction": 0.0
        },
        # Refractory Asthenosphere
        13: {
            "mat_name": mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type": types.ultramafic_mantle_refactory,
            "mat_domain": domains.general_domain,
            "red_fraction": 0.0,
            "green_fraction": 255/255,
            "blue_fraction": 0.0
        },
        # Solidified Gabbro
        14: {
            "mat_name": mat_names.oceanic_gabbroic_crust,
            "mat_type": types.solidified_gabbro,
            "mat_domain": domains.general_domain,
            "red_fraction": 227/255,
            "green_fraction": 227/255,
            "blue_fraction": 227/255
        },
        # Partially Molten Gabbro
        15: {
            "mat_name": mat_names.oceanic_gabbroic_crust,
            "mat_type": types.solidified_gabbro_partially_molten,
            "mat_domain": domains.general_domain,
            "red_fraction": 255/255,
            "green_fraction": 125/255,
            "blue_fraction": 255/255
        },
        # Extracted Gabbroic Magma
        16: {
            "mat_name": mat_names.gabbroic_magma,
            "mat_type": types.extracted_gabbroic_magma,
            "mat_domain": domains.general_domain,
            "red_fraction": 255/255,
            "green_fraction": 100/255,
            "blue_fraction": 100/255
        },
        # Layered Solidified Gabbroic Crust
        17: {
            "mat_name": mat_names.layered_gabbroic_crust,
            "mat_type": types.solidified_layered_gabbro,
            "mat_domain": domains.general_domain,
            "red_fraction": 175/255,
            "green_fraction": 175/255,
            "blue_fraction": 175/255
        },
        # Layered Partially Molten Gabbroic Crust
        18: {
            "mat_name": mat_names.layered_gabbroic_crust,
            "mat_type": types.solidified_layered_gabbro_partially_molten,
            "mat_domain": domains.general_domain,
            "red_fraction": 255/255,
            "green_fraction": 0/255,
            "blue_fraction": 255/255
        },
        # Layered Gabbroic Magma
        19: {
            "mat_name": mat_names.layered_gabbroic_magma,
            "mat_type": types.extracted_layered_gabbroic_magma,
            "mat_domain": domains.general_domain,
            "red_fraction": 255/255,
            "green_fraction": 0/255,
            "blue_fraction": 0/255
        },
        # Extruded Gabbroic Magma (Lava)
        20: {
            "mat_name": mat_names.gabbroic_magma,
            "mat_type": types.extruded_gabbroic_magma,
            "mat_domain": domains.general_domain,
            "red_fraction": 255/255,
            "green_fraction": 255/255,
            "blue_fraction": 0/255
        },
        # Solidified Extruded Gabbroic Magma (Basalt)
        21: {
            "mat_name": mat_names.oceanic_gabbroic_crust,
            "mat_type": types.solidified_basalt,
            "mat_domain": domains.general_domain,
            "red_fraction": 168/255,
            "green_fraction": 45/255,
            "blue_fraction": 0/255
        }
    }

    return materials_input_dict, lib.path
