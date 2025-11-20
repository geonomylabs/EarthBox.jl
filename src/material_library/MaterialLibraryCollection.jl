module MaterialLibraryCollection

include("utils/MaterialCollectionManager.jl")
include("registries/__init__.jl")

using EarthBox
import YAML
import DataStructures: OrderedDict
import .MaterialCollectionManager: MaterialCollection
using .BenchmarksRegistry
using .LithosphericDeformationRegistry: LithosphericDeformationCollections
using .Gerya2019Registry
using .ViscoElastoPlasticRegistry
using .SimpleMantleMeltingRegistry
using .SeafloorSpreadingRegistry
using .SandboxRegistry

function get_material_names_list_string(collection::MaterialCollection)::String
    material_names = collection.materials
    return join(["- `$(material_name)`" for material_name in material_names], "\n")
end

function make_material_table(
    collection::MaterialCollection;
    target_collection_name::String,
    material_name::String
)::String
    collection_path = collection.path
    yaml_data = YAML.load_file(collection_path, dicttype=OrderedDict)
    collection_key = nothing
    for key in keys(yaml_data)
        if key == target_collection_name
            collection_key = key
            break
        end
    end
    if collection_key === nothing
        error("Material collection '$target_collection_name' not found in YAML file")
    end
    materials_data = yaml_data[collection_key]
    table_lines = String[]
    push!(table_lines, "| Property Name | Value | Units | Description |")
    push!(table_lines, "|:--------------|:------|:------|:------------|")
    if !haskey(materials_data, material_name)
        @warn "Material $material_name' not found in YAML file"
        return "No table available for material $material_name"
    end
    material_properties = materials_data[material_name]
    for (prop_name, prop_data) in material_properties
        if isa(prop_data, Vector) && length(prop_data) >= 3
            value = prop_data[1]
            units = prop_data[2]
            description = prop_data[3]
            push!(table_lines, "| `$(prop_name)` | `$(value)` | `$(units)` | `$(description)` |")
        end
    end
    return join(table_lines, "\n")
end

function make_material_collection_table(
    collection::MaterialCollection;
    target_collection_name::String
)::Tuple{String, String}
    material_names = collection.materials
    collection_path = collection.path
    yaml_data = YAML.load_file(collection_path, dicttype=OrderedDict)
    collection_key = nothing
    for key in keys(yaml_data)
        if key == target_collection_name
            collection_key = key
            break
        end
    end
    main_description = ""
    for key in keys(yaml_data)
        if key == "description"
            main_description = yaml_data[key]
        end
    end
    if collection_key === nothing
        error("Material collection '$target_collection_name' not found in YAML file")
    end
    materials_data = yaml_data[collection_key]
    table_lines = String[]
    # Add table header
    push!(table_lines, "| Material Name | Property Name | Value | Units | Description |")
    push!(table_lines, "|--------------:|:--------------|:------|:------|:------------|")
    for material_name in material_names
        if !haskey(materials_data, material_name)
            @warn "Material $material_name' not found in YAML file"
            continue
        end
        push!(table_lines, "| **`$(material_name)`** | | | | |")
        material_properties = materials_data[material_name]
        for (prop_name, prop_data) in material_properties
            if isa(prop_data, Vector) && length(prop_data) >= 3
                value = prop_data[1]
                units = prop_data[2]
                description = prop_data[3]
                push!(table_lines, "| | `$(prop_name)` | `$(value)` | `$(units)` | `$(description)` |")
            end
        end
    end
    return join(table_lines, "\n"), main_description
end

"""
    MaterialLibrary

The material library struct is composed of material collections that define 
material names associated with a material collection and the path to a yamal 
formatted file with material properties associated with each material name in 
the collection.

# Fields
- `lithospheric_deformation`::[`LithosphericDeformationCollections`](@ref):
    - Collection of material collections for lithospheric deformation models.
    - See yamal files in src/material_library/registries/lithospheric_deformation.
- `sandbox`::[`MaterialCollection`](@ref):
    - Material names and material collection file path for sandbox experiments.
    - See yamal files in src/material_library/registries/sandbox.
- `benchmarks`::[`MaterialCollection`](@ref):
    - Material names and material collection file path for benchmark models.
    - See yamal files in src/material_library/registries/benchmarks.
- `gerya2019`::[`MaterialCollection`](@ref):
    - Material names and material collection file path for the models of Gerya (2019).
    - See yamal files in src/material_library/registries/gerya2019.
- `viscoelastoplastic`::[`MaterialCollection`](@ref):
    - Material names and material collection file path for viscoelasto-plastic models.
    - See yamal files in src/material_library/registries/viscoelastoplastic.
- `simple_mantle_melting`::[`MaterialCollection`](@ref):
    - Material names and material collection file path for the simple mantle melting model.
    - See yamal files in src/material_library/registries/simple_mantle_melting.
- `seafloor_spreading`::[`MaterialCollection`](@ref):
    - Material names and material collection file path for the seafloor spreading model.
    - See yamal files in src/material_library/registries/seafloor_spreading.
"""
struct MaterialLibrary
    lithospheric_deformation::LithosphericDeformationCollections
    sandbox::MaterialCollection
    benchmarks::MaterialCollection
    gerya2019::MaterialCollection
    viscoelastoplastic::MaterialCollection
    simple_mantle_melting::MaterialCollection
    seafloor_spreading::MaterialCollection
end

function MaterialLibrary()
    return MaterialLibrary(
        LithosphericDeformationCollections(),
        SandboxRegistry.get_material_collection(),
        BenchmarksRegistry.get_material_collection(),
        Gerya2019Registry.get_material_collection(),
        ViscoElastoPlasticRegistry.get_material_collection(),
        SimpleMantleMeltingRegistry.get_material_collection(),
        SeafloorSpreadingRegistry.get_material_collection()
    )
end

end # module 