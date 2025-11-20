# Material Input

Material inputs define information associated with each material used in the model
and must be provided as an argument to [MarkerMaterials.initialize!](@ref Markers.MarkerMaterials.initialize!)
either in the form of a path to a yamal formatted file or Julia dictionary.

Material inputs must include the following information for each material:
- Integer ID greater than or equal to 1
- Material name (`mat_name`) from the user selected material collection library file. See 
   [Material Collection Files](@ref) and [MaterialLibrary](@ref) for information on how to select
   and query a material collection library file packaged with EarthBox or build one from scratch.
- Material type name (`mat_type`) that defines evolving genetic characteristics. See 
   [Material Types](@ref) for more information about material type.
- Material domain name (`mat_domain`) that defines the physical domain associated with the material.
   See [Material Domains](@ref) for more information about material domains.
- red-green-blue colors used in EarthBox plotting functions.

# Material Input Files

Material input files are .yaml files with the following structure:

```
name: user_defined_name_of_material_model
description: 'Add text to describe this material model'
1:
  mat_name: [MaterialCollectionName1, None, ' Material name from material library ']
  mat_type: [MaterialTypeName1, None, ' Material type ']
  mat_domain: [MaterialDomainName1, None, ' Material domain ']
  red_fraction: [1.0, None, ' red_fraction']
  green_fraction: [1.0, None, ' green_fraction']
  blue_fraction: [1.0, None, ' blue_fraction ']
2:
  mat_name: [MaterialCollectionName2, None, ' Material name from material library ']
  mat_type: [MaterialTypename2, None, ' Material type ']
  mat_domain: [MaterialDomainName2, None, ' Material domain ']
  red_fraction: [0.0, None, ' red_fraction']
  green_fraction: [1.0, None, ' green_fraction']
  blue_fraction: [1.0, None, ' blue_fraction ']
3: ...
```

where `MaterialCollectionName1` and `MaterialCollection2` are material names from the selected material 
collection library file ([Material Collection Files](@ref)), `MaterialTypename1` and `MaterialTypeName2`
are genetic material type names ([Material Types](@ref)), `MaterialDomainName1` and `MaterialDomainname2`
are material domain names ([Material Domains](@ref)), and `red_fraction`, `green_fraction` and
`blue_fraction` are RGB color values. The information provided for each information type is in the form
of a list with entries for value, units and description. The units are `None` for the current 
implementation and the description can be used for general documentation. 

!!! warning "material ID's"
    Integer ID's must be greater than or equal to 1.

!!! info "accessing valid names"
    Valid material, type and domain names and definitions can be accessed using the instructions from
    [MaterialLibrary](@ref), [Material Types](@ref) and [Material Domains](@ref).


# Material Input Dictionaries

Material input can also be provided as a Julia dictionary. This gives the user the ability to
programmatically access valid material, genetic type names and domain names. The following example
uses a material collection library file from EarthBox to build a material input dictionary for a 
model with a sticky air/water layer, continental lithosphere, lateral lithospheric strong zones and 
melting processes including melt extraction and melt extrusion.

```julia
using EarthBox

MaterialDictType = EarthBoxDtypes.MaterialDictType;
MaterialsDictType = EarthBoxDtypes.MaterialsDictType;

collection = MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1;
mat_names = collection.materials;
types = MaterialTypesRegistry();
domains = materialDomainsRegistry()

material_input_dict = MaterialsDictType(
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
        "red_fraction" => 0.89411765,
        "green_fraction" => 0.58039216,
        "blue_fraction" => 0.28627451,
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
        "mat_domain" => domains.upper_mantle_lithosphere,
        "red_fraction" => 0.0,
        "green_fraction" => 153/255,
        "blue_fraction" => 153/255,
    ),
    Int16(9) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_continental_lithosphere,
        "mat_type" => types.ultramafic_mantle_fertile,
        "mat_domain" => domains.middle_mantle_lithosphere,
        "red_fraction" => 0.0,
        "green_fraction" => 153/255,
        "blue_fraction" => 153/255,
    ),
    Int16(10) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_continental_lithosphere,
        "mat_type" => types.ultramafic_mantle_fertile,
        "mat_domain" => domains.lower_mantle_lithosphere,
        "red_fraction" => 0.0,
        "green_fraction" => 153/255,
        "blue_fraction" => 153/255,
    ),
    Int16(11) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_continental_lithosphere_strong_zone,
        "mat_type" => types.ultramafic_mantle_fertile,
        "mat_domain" => domains.lithospheric_mantle_strong_zone,
        "red_fraction" => 0.0,
        "green_fraction" => 163/255,
        "blue_fraction" => 163/255,
    ),
    # Asthenosphere
    Int16(12) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
        "mat_type" => types.ultramafic_mantle_fertile,
        "mat_domain" => domains.asthenosphere,
        "red_fraction" => 0.0,
        "green_fraction" => 200/255,
        "blue_fraction" => 153/255,
    ),
    # Partially Molten Asthenosphere
    Int16(13) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
        "mat_type" => types.ultramafic_mantle_partially_molten,
        "mat_domain" => domains.general_domain,
        "red_fraction" => 255/255,
        "green_fraction" => 150/255,
        "blue_fraction" => 0.0,
    ),
    # Refractory Asthenosphere
    Int16(14) => MaterialDictType(
        "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
        "mat_type" => types.ultramafic_mantle_refactory,
        "mat_domain" => domains.general_domain,
        "red_fraction" => 0.0,
        "green_fraction" => 255/255,
        "blue_fraction" => 0.0,
    ),
    # Solidified Gabbro
    Int16(15) => MaterialDictType(
        "mat_name" => mat_names.oceanic_gabbroic_crust,
        "mat_type" => types.solidified_gabbro,
        "mat_domain" => domains.general_domain,
        "red_fraction" => 227/255,
        "green_fraction" => 227/255,
        "blue_fraction" => 227/255,
    ),
    # Partially Molten Gabbro
    Int16(16) => MaterialDictType(
        "mat_name" => mat_names.oceanic_gabbroic_crust,
        "mat_type" => types.solidified_gabbro_partially_molten,
        "mat_domain" => domains.general_domain,
        "red_fraction" => 255/255,
        "green_fraction" => 125/255,
        "blue_fraction" => 255/255,
    ),
    # Extracted Gabbroic Magma
    Int16(17) => MaterialDictType(
        "mat_name" => mat_names.gabbroic_magma,
        "mat_type" => types.extracted_gabbroic_magma,
        "mat_domain" => domains.general_domain,
        "red_fraction" => 255/255,
        "green_fraction" => 100/255,
        "blue_fraction" => 100/255,
    ),
    # Layered Solidified Gabbroic Crust
    Int16(18) => MaterialDictType(
        "mat_name" => mat_names.layered_gabbroic_crust,
        "mat_type" => types.solidified_layered_gabbro,
        "mat_domain" => domains.general_domain,
        "red_fraction" => 175/255,
        "green_fraction" => 175/255,
        "blue_fraction" => 175/255,
    ),
    # Layered Partially Molten Gabbroic Crust
    Int16(19) => MaterialDictType(
        "mat_name" => mat_names.layered_gabbroic_crust,
        "mat_type" => types.solidified_layered_gabbro_partially_molten,
        "mat_domain" => domains.general_domain,
        "red_fraction" => 255/255,
        "green_fraction" => 0/255,
        "blue_fraction" => 255/255,
    ),
    # Layered Gabbroic Magma
    Int16(20) => MaterialDictType(
        "mat_name" => mat_names.layered_gabbroic_magma,
        "mat_type" => types.extracted_layered_gabbroic_magma,
        "mat_domain" => domains.general_domain,
        "red_fraction" => 255/255,
        "green_fraction" => 0/255,
        "blue_fraction" => 0/255,
    ),
    # Extruded Gabbroic Magma (Lava)
    Int16(21) => MaterialDictType(
        "mat_name" => mat_names.gabbroic_magma,
        "mat_type" => types.extruded_gabbroic_magma,
        "mat_domain" => domains.general_domain,
        "red_fraction" => 255/255,
        "green_fraction" => 255/255,
        "blue_fraction" => 0/255,
    ),
    # Solidified Extruded Gabbroic Magma (Basalt)
    Int16(22) => MaterialDictType(
        "mat_name" => mat_names.oceanic_gabbroic_crust,
        "mat_type" => types.solidified_basalt,
        "mat_domain" => domains.general_domain,
        "red_fraction" => 168/255,
        "green_fraction" => 45/255,
        "blue_fraction" => 0/255,
    )
)
```

!!! note "integer ID's"
    Integer ID's must be defined with Int16 to maintain compatibility with material
    ID marker arrays.

