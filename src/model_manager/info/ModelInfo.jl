module ModelInfo

using YAML
using DataStructures: OrderedDict
import EarthBox.ModelDataContainer: ModelData
import Statistics

""" Print model array or parameter information.

Inputs
------
`model:: ModelData`: The model data container
`output_dir_path:: String`: Directory path for output
`collection_type:: String`: Either "arrays" or "parameters"
"""
function export_collection_information(
    model::ModelData,
    output_dir_path::String,
    collection_type::String
)::Nothing
    collection_dict = get_collection_dict(model, collection_type)
    group_object_list, group_object_paths_list = get_group_object_list(collection_dict)
    
    info_dict = make_info_dict(
        collection_type,
        group_object_list,
        group_object_paths_list
    )
    sorted_info_dict = OrderedDict{String, Any}(k => info_dict[k] for k in sort(collect(keys(info_dict))))
    save_information_to_file(output_dir_path, sorted_info_dict, collection_type)
    return nothing
end

""" Build a list of information for each collection.

A collection is a collection of either arrays or parameters associated with a model data 
container that can be accessed within the code via model.container_name.collection_type 
where container_name is the name of the main model data container.

Dictionary keys are defined as a string defining the full reference to the array collection:
    model.container_name.arrays.

Each dictionary entry is a list where the first entry is equal to the instance of the array 
collection class and the second entry is a list of array group names:
    [arrays_collection_obj, [group_name1, group_name2,....]]

If there are no group names then this list will be empty.
"""
function get_collection_dict(
    model::ModelData,
    collection_type::Union{String, Symbol}
)::OrderedDict{String, Vector{Any}}
    collection_dict = OrderedDict{String, Vector{Any}}()
    
    if collection_type ∉ ["arrays", "parameters", :arrays, :parameters]
        error("Collection type must be arrays or parameters or :arrays or :parameters")
    end
    
    for container_name in propertynames(model)
        container_obj = getproperty(model, container_name)
        container_property_names = propertynames(container_obj)
        if Symbol(collection_type) in container_property_names
            collection_obj = getproperty(container_obj, Symbol(collection_type))
            collection_name_str = "model.$(container_name).$(collection_type)"
            group_name_list = get_group_names(collection_obj, collection_type)
            if length(group_name_list) > 0
                collection_dict[collection_name_str] = [collection_obj, group_name_list]
            else
                collection_dict["model.$(container_name)"] = [container_obj, String[collection_type]]
            end
        end
    end
    return collection_dict
end

""" Get group names associated with the collection object.
"""
function get_group_names(collection_obj::Any, collection_type::String)::Vector{String}
    group_name_list = String[]
    
    for group_name in propertynames(collection_obj)
        grp_obj = getproperty(collection_obj, group_name)
        if collection_type == "arrays"
            isarray = isarrayobject(grp_obj)
            if !isarray && group_name != :obj_list
                push!(group_name_list, String(group_name))
            end
        elseif collection_type == "parameters"
            isparameter = isparameterobject(grp_obj)
            if !isparameter && group_name != :obj_list
                push!(group_name_list, String(group_name))
            end
        end
    end
    return group_name_list
end

""" Check if the object is an array object.
"""
function isarrayobject(obj::Any)::Bool
    return hasproperty(obj, :array)
end

""" Check if the object is a parameter object.
"""
function isparameterobject(obj::Any)::Bool
    return hasproperty(obj, :value)
end

""" Make a list of group instances for the collection.
"""
function get_group_object_list(
    collection_dict::OrderedDict{String, Vector{Any}}
)::Tuple{Vector{Any}, Vector{String}}
    group_object_list = Any[]
    group_object_paths_list = String[]
    
    for (collection_name, collection_list) in collection_dict
        collection_obj = collection_list[1]
        group_name_list = collection_list[2]
        
        if !isempty(group_name_list)
            for group_name in group_name_list
                group_obj = getproperty(collection_obj, Symbol(group_name))
                push!(group_object_list, group_obj)
                push!(group_object_paths_list, "$(collection_name).$(group_name)")
            end
        else
            push!(group_object_list, collection_obj)
            push!(group_object_paths_list, collection_name)
        end
    end
    return group_object_list, group_object_paths_list
end

""" Make a dictionary of target objects (arrays or parameters).
"""
function make_info_dict(
    collection_type::String,
    group_object_list::Vector{Any},
    group_object_paths_list::Vector{String}
)::OrderedDict{String, Any}
    info_dict = OrderedDict{String, Any}()
    # ignore attribute names in group object that are not associated with array
    # or parameter objects.
    names_to_ignore = ["group_name"]
    
    for (i, group_obj) in enumerate(group_object_list)
        group_obj_path = group_object_paths_list[i]
        for attr_name in propertynames(group_obj)
            target_name = String(attr_name)
            if !(target_name in ["obj_list"]) && !(target_name in names_to_ignore)
                target_obj = getproperty(group_obj, attr_name)
                collect_info!(
                    collection_type,
                    target_name,
                    target_obj,
                    group_obj_path,
                    info_dict
                )
            end
        end
    end
    return info_dict
end

""" Collect array or parameter information and add to info_dict.
"""
function collect_info!(
    collection_type::String,
    target_name::String,
    target_obj::Any,
    group_obj_path::String,
    info_dict::OrderedDict{String, Any}
)::Nothing
    name = target_name
    
    units = hasproperty(target_obj, :units) ? target_obj.units : "NA"
    description = hasproperty(target_obj, :description) ? target_obj.description : "None"
    
    array_shape = hasproperty(target_obj, :array) ? size(target_obj.array) : nothing

    obj_min = "NA"
    obj_max = "NA"
    obj_value = "NA"
    obj_mean = "NA"
    if hasproperty(target_obj, :array)
        obj_type = typeof(target_obj.array)
        obj_min = minimum(target_obj.array)
        obj_max = maximum(target_obj.array)
        obj_mean = Statistics.mean(target_obj.array)
    elseif hasproperty(target_obj, :value)
        obj_type = typeof(target_obj.value)
        obj_value = target_obj.value
    else
        obj_type = nothing
        obj_value = nothing
    end
    
    units_list = isa(units, String) ? [units] : collect(units)
    description_list = isa(description, String) ? [description] : collect(description)
    
    if length(units_list) == 1
        units = units_list[1]
    else
        units = units_list
    end

    if length(description_list) == 1
        description = description_list[1]
    else
        description = description_list
    end

    if collection_type == "arrays"
        if length(array_shape) == 1
            shape_str = string(array_shape[1])
        else
            shape_str = string(array_shape[1], ",", array_shape[2])
        end
        info_dict[name] = OrderedDict{String, Any}(
            "name" => name,
            "reference" => "$(group_obj_path).$(name).array",
            "units" => units,
            "shape" => shape_str,
            "type" => string(obj_type),
            "description" => description,
            "min" => obj_min,
            "max" => obj_max,
            "mean" => obj_mean
        )
    else
        info_dict[name] = OrderedDict{String, Any}(
            "name" => name,
            "reference" => "$(group_obj_path).$(name).value",
            "value" => obj_value,
            "units" => units,
            "type" => string(obj_type),
            "description" => description
        )
    end
    return nothing
end

""" Write array info to output file.
"""
function save_information_to_file(
    output_dir_path::String,
    info_dict::OrderedDict{String, Any},
    collection_type::String
)::Nothing
    file_path = joinpath(output_dir_path, "$(collection_type).yml")
    make_yaml_file(file_path, info_dict)
    return nothing
end

""" Write dictionary to yaml file.
"""
function make_yaml_file(yaml_file_path::String, data_dict::OrderedDict{String, Any})::Nothing
    open(yaml_file_path, "w") do io
        YAML.write(io, data_dict)
    end
    return nothing
end

end # module 