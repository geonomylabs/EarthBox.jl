module MaterialsContainer

include("properties/MaterialProperties.jl")
include("material/MaterialContainer.jl")
include("loader/LoadMaterialArrays.jl")
include("loader/LoadMaterialDicts.jl")
include("reader/MaterialReader.jl")
include("utils/MaterialsStateContainer.jl")

import CairoMakie
import EarthBox.ModelDataContainer: ModelData
import EarthBox.MaterialColorsContainer: MaterialColors
import EarthBox.MaterialColorsContainer: get_colorbar_ticks_for_material_colors
import EarthBox.PrintFuncs: print_info
import EarthBox.Parameters: ParameterInt, ParameterFloat, ParameterStr
import EarthBox.EarthBoxDtypes: MaterialsDictType
import EarthBox.OptionTools: OptionState
import .MaterialProperties
import .MaterialContainer: Material
import .MaterialContainer: check_and_clean_materials_input_dict
import .LoadMaterialArrays: load_material_arrays
import .LoadMaterialArrays: MaterialCollectionDictType
import .LoadMaterialDicts: update_material_id_dicts
import .MaterialReader: read_materials_library
import .MaterialReader: read_materials_input
import .MaterialReader: build_materials_dict
import .MaterialsStateContainer: MaterialsState
import .MaterialsStateContainer: update

""" Materials struct.

# Attributes
- material_library_file_path::Union{String, Nothing}
    - Path to the material library file containing material properties for
    each material name.

- material_model_file_path::Union{String, Nothing}
    - Path to the material model file containing references to the material
    names stored in the material library for each material id in the
    model. This file also contains the material type, domain and color
    for each material id in the model.

- materials::Dict{String, Material}
    - Dictionary of material objects where the key in material id and the
    value is an instance of a material class. The Material class has
    attributes for different material property class instances that
    provide access to the property values loaded from the materials_dict.

- materials_input_dict::Union{MaterialsDictType, Nothing}
    - Dictionary of material input dictionaries. Each main key is a material
    id and the value is a dictionary of material input properties including
    the material name, type, domain and red, green and blue color values.

- materials_dict::MaterialsDictType
    - Dictionary of material dictionaries. Each main key is a material id
    used in the model and the value is a dictionary of material parameters
    where keys are material parameter names and values are parameter values.

    For example:
        materials_dict[matid][parameter_name] = parameter_value

    This dictionary includes the entries from the materials_input_dict
    including name, type, domain and red, green and blue color values.

- materials_library_dict::MaterialsDictType
    - Dictionary of materials from the material library. Each main key is a
    material name and the value is a dictionary of material properties.

    For example:
        materials_library_dict[material_name][parameter_name] = parameter_value

- itype_mat::Union{Int64, Nothing}
    - Material type integer id.

- material_option::Union{OptionState, Nothing}
    - Material model option.

- colors::MaterialColors
    - Struct with material color information.

"""
mutable struct Materials
    material_library_file_path::Union{String, Nothing}
    material_model_file_path::Union{String, Nothing}
    materials::MaterialCollectionDictType
    materials_input_dict::Union{MaterialsDictType, Nothing}
    materials_dict::MaterialsDictType
    materials_library_dict::MaterialsDictType
    itype_mat::Union{Int64, Nothing}
    material_option::Union{OptionState, Nothing}
    colors::MaterialColors

    function Materials()
        new(
            nothing,
            nothing,
            MaterialCollectionDictType(),
            nothing,
            MaterialsDictType(),
            MaterialsDictType(),
            nothing,
            nothing,
            MaterialColors()
        )
    end
end

function set_material_option!(
    materials::Materials,
    itype_mat::Int64,
    material_option::OptionState
)::Nothing
    materials.itype_mat = itype_mat
    materials.material_option = material_option
    return nothing
end

function load_materials_after_checks!(
    materials::Materials,
    material_library_file_path::Union{String, Nothing},
    material_model_file_path::Union{String, Nothing},
    materials_input_dict::Union{MaterialsDictType, Nothing}=nothing,
    materials_state::Union{MaterialsState, Nothing}=nothing
)::Union{MaterialsState, Nothing}
    materials_input_exists = is_material_input(material_model_file_path, materials_input_dict)
    if material_library_file_path !== nothing && materials_input_exists
        load!(
            materials,
            material_library_file_path=material_library_file_path,
            material_model_file_path=material_model_file_path,
            materials_input_dict=materials_input_dict
        )
        if materials_state !== nothing
            materials_state =update(materials_state, "Loaded")
        end
    end
    return materials_state
end

function is_material_input(
    material_model_file_path::Union{String, Nothing},
    materials_input_dict::Union{MaterialsDictType, Nothing}
)::Bool
    materials_input_exists = true
    if material_model_file_path === nothing && materials_input_dict === nothing
        materials_input_exists = false
    end
    return materials_input_exists
end

function load!(
    materials::Materials;
    model::Union{ModelData, Nothing}=nothing,
    material_library_file_path::Union{String, Nothing}=nothing,
    material_model_file_path::Union{String, Nothing}=nothing,
    materials_input_dict::Union{MaterialsDictType, Nothing}=nothing
)::Nothing
    materials.material_library_file_path = material_library_file_path
    materials.material_model_file_path = material_model_file_path
    materials.materials_input_dict = materials_input_dict
    print_paths(materials)
    
    set_materials_library_dict!(materials)
    set_materials_input_dict!(materials)
    materials.materials_dict, _ = build_materials_dict(
        materials.materials_library_dict,
        materials.materials_input_dict
    )
    set_materials!(materials)
    # This chunk is not used for post-processing cases where model is not provided
    if length(materials.materials) > 0 && model !== nothing
        update_material_arrays_and_dicts!(materials, model)
    end
    return nothing
end

function print_paths(materials::Materials)::Nothing
    print_info("Material library file path: $(materials.material_library_file_path)", level=1)
    print_info("Material model file path: $(materials.material_model_file_path)", level=1)
    return nothing
end

function set_materials_library_dict!(materials::Materials)::Nothing
    if materials.material_library_file_path === nothing
        error("The path to the material library must be defined.")
    end
    materials.materials_library_dict, _ = read_materials_library(
        materials.material_library_file_path
    )
    return nothing
end

function set_materials_input_dict!(materials::Materials)::Nothing
    if materials.material_model_file_path !== nothing && materials.material_model_file_path != ""
        materials.materials_input_dict, _ = read_materials_input(
            materials.material_model_file_path
        )
    elseif materials.materials_input_dict !== nothing
        materials.materials_input_dict = check_and_clean_materials_input_dict(
            materials.materials_input_dict,
            materials.material_option
        )
    else
        error("A path or dictionary for material inputs were not provided. Adjust inputs.")
    end
    return nothing
end

function set_materials!(materials::Materials)::Nothing
    materials_dict = Dict{Union{String, Int16}, Material}()
    for (matid, material_dict) in materials.materials_dict
        if typeof(matid) == String
            matid = parse(Int16, matid)  # Convert string to Int16
        else
            matid = convert(Int16, matid)  # Convert Int64 to Int16
        end
        materials_dict[matid] = Material(material_dict)
    end
    materials.materials = materials_dict
    return nothing
end

function update_material_arrays_and_dicts!(
    materials::Materials,
    model::ModelData
)::Nothing
    load_material_arrays(materials.materials, model)
    update_material_id_dicts(materials.materials, model)
    update_material_colors!(materials, model)
    return nothing
end

function print_materials_input_dict(materials::Materials)::Nothing
    dict = materials.materials_input_dict
    print_info("", level=1)
    print_info("Materials input dict:", level=1)
    print_info("", level=1)
    for (key, inner_dict) in dict
        print_info("Key: $key", level=2)
        for (inner_key, value) in inner_dict
            print_info("  $inner_key: $value", level=3)
        end
        print_info("", level=1)
    end
    return nothing
end

function print_materials_dict(materials::Materials)::Nothing
    dict = materials.materials_dict
    print_info("", level=1)
    print_info("Materials dict:", level=1)
    print_info("", level=1)
    for (key, inner_dict) in dict
        print_info("Key: $key", level=2)
        for (inner_key, value) in inner_dict
            print_info("  $inner_key: $value", level=3)
        end
        print_info("", level=1)
    end
    return nothing
end

function print_materials_library_dict(materials::Materials)::Nothing
    dict = materials.materials_library_dict
    print_info("", level=1)
    print_info("Materials library dict:", level=1)
    print_info("", level=1)
    for (key, inner_dict) in dict
        print_info("Key: $key", level=2)
        for (inner_key, value) in inner_dict
            print_info("  $inner_key: $value", level=3)
        end
        print_info("", level=1)
    end
    return nothing
end

function print_materials(materials::Materials)::Nothing
    open("materials.dat", "w") do io
        dict = materials.materials
        print_info("", level=1)
        print_info("Materials dict:", level=1)
        print_info("", level=1)
        for (key, material) in dict
            print_info("Key: $key", level=2)
            for field in fieldnames(typeof(material))
                value = getfield(material, field)
                print_info("$field: $value", level=3)
            end
            print_info("", level=1)
        end
    end
    return nothing
end

function update_material_colors!(materials::Materials, model::ModelData)::Nothing
    materials.colors.color_map, materials.colors.n_bin = make_color_map_comp(materials)
    materials.colors.labels = make_custom_labels(materials)
    materials.colors.ticks = get_colorbar_ticks_for_material_colors(materials.colors.n_bin)
    materials.colors.colorrange = (1.0 - 0.5, Float64(materials.colors.n_bin) + 0.5)
    # Now update the persistent model data container
    model.materials.colors = materials.colors
    return nothing
end

function make_color_map_comp(materials::Materials)::Tuple{Any, Int64}
    materials_dict = materials.materials
    material_ids = collect(keys(materials_dict))
    ncolors = maximum(material_ids)
    colors = Vector{CairoMakie.RGB{Float64}}()
    
    for matid in 1:ncolors
        if matid in material_ids
            rgb = materials_dict[matid].rgb
            red_fraction = rgb.red_fraction.value
            green_fraction = rgb.green_fraction.value
            blue_fraction = rgb.blue_fraction.value
        else
            red_fraction = 0.0
            green_fraction = 0.0
            blue_fraction = 0.0
        end
        color = CairoMakie.RGB(red_fraction, green_fraction, blue_fraction)
        push!(colors, color)
    end
    
    n_bin = ncolors
    if length(colors) < 2
        push!(colors, CairoMakie.RGB(1.0, 1.0, 1.0))
        n_bin += 1
    end
    
    cm = CairoMakie.cgrad(colors; categorical = true)
    return cm, n_bin
end

function make_custom_labels(materials::Materials)::Vector{String}
    materials_dict = materials.materials_dict
    # Create sorted list of material IDs by converting to Int16 and then sorting
    # But first check if the keys are strings or integers
    dict_keys = collect(keys(materials_dict))
    matid_test = dict_keys[1]
    if isa(matid_test, String)
        matids = sort([parse(Int16, matid) for matid in keys(materials_dict)])
    else
        matids = sort([matid for matid in keys(materials_dict)])
    end
    custom_labels = String[]
    for matid in matids
        # Convert matid back to a string since this is the expected type
        if isa(matid_test, String)
            matid = string(matid)
        end
        params = materials_dict[matid]
        mat_name = params["mat_name"]
        mat_name = replace_underscores_with_spaces(mat_name)
        mat_type = params["mat_type"]
        mat_type = add_spaces_in_front_of_capitals(mat_type)
        mat_domain = params["mat_domain"]
        mat_domain = add_spaces_in_front_of_capitals(mat_domain)
        spacer = get_spacer(matid)
        label = "$matid: $mat_name \n$spacer($mat_type : $mat_domain)"
        push!(custom_labels, label)
    end
    return custom_labels
end

function get_spacer(matid::Union{String, Int16})::String
    matid_str = string(matid)
    if length(matid_str) == 1
        spacer = "    "  # 4 spaces
    elseif length(matid_str) == 2
        spacer = "     "  # 5 spaces
    else
        spacer = "      "  # 6 spaces
    end
    return spacer
end

function add_spaces_in_front_of_capitals(string::String)::String
    result = ""
    for (i, char) in enumerate(string)
        if i > 1 && isuppercase(char)
            result *= " "
        end
        result *= char
    end
    return result
end

function replace_underscores_with_spaces(string::String)::String
    return replace(string, "_" => " ")
end

end # module