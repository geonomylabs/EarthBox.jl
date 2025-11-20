module ObjDictUtils

import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D, AbstractEarthBoxArray2D
import EarthBox.EarthBoxDtypes: AbstractParameter
import EarthBox.EarthBoxDtypes: AbstractParameterGroup, AbstractArrayGroup
import EarthBox.EarthBoxDtypes: AbstractParameterCollection, AbstractArrayCollection
import EarthBox.EarthBoxDtypes: CollectionContainer
import EarthBox.EarthBoxDtypes: ObjDictType

""" 
Adds all collections from a CollectionContainer to the provided dictionary. 
Iterates through all fields of the CollectionContainer and adds any collections 
that are either AbstractParameterCollection or AbstractArrayCollection.

# Arguments
- obj_dict::Dict{String, T}: Dictionary to add collections to
- collection_container::CollectionContainer: Container holding collections to add

# Returns
- Nothing
"""
function add_collections_to_dict!(
    obj_dict::Dict{String, T},
    collection_container::CollectionContainer
)::Nothing where T
    for name in fieldnames(typeof(collection_container))
        collection = getfield(collection_container, name)
        if (
            collection isa AbstractParameterCollection || 
            collection isa AbstractArrayCollection
        )
            add_collection_to_dict!(obj_dict, collection)
        elseif check_for_earthbox_object(collection)
            manage_check_and_add(collection, obj_dict)
        end
    end
    return nothing
end

""" 
Adds all groups from a collection to the provided dictionary. Iterates through 
all fields of the collection and adds any groups that are either 
AbstractParameterGroup or AbstractArrayGroup.

# Arguments
- obj_dict::Dict{String, T}: Dictionary to add groups to
- collection::Union{AbstractParameterCollection, AbstractArrayCollection}: 
    Collection containing groups to add

# Returns
- Nothing
"""
function add_collection_to_dict!(
    obj_dict::Dict{String, T},
    collection::Union{AbstractParameterCollection, AbstractArrayCollection}
)::Nothing where T
    for name in fieldnames(typeof(collection))
        group = getfield(collection, name)
        if (
            group isa AbstractParameterGroup || 
            group isa AbstractArrayGroup
        )
            add_objects_to_dict!(obj_dict, group)
        elseif check_for_earthbox_object(group)
            manage_check_and_add(group, obj_dict)
        end
    end
    return nothing
end

"""
Adds all EarthBox objects from a group to the provided dictionary. Iterates 
through all fields of the group and adds any objects that are AbstractParameter, 
AbstractEarthBoxArray1D, or AbstractEarthBoxArray2D and have a name property.

# Arguments
- obj_dict::Dict{String, T}: Dictionary to add objects to
- group::Union{AbstractParameterGroup, AbstractArrayGroup}: Group containing 
    objects to add

# Returns
- Nothing
"""
function add_objects_to_dict!(
    obj_dict::Dict{String, T},
    group::Union{AbstractParameterGroup, AbstractArrayGroup}
)::Nothing where T
    for name in fieldnames(typeof(group))
        earthbox_object = getfield(group, name)
        manage_check_and_add(earthbox_object, obj_dict)
    end
    return nothing
end

function manage_check_and_add(
    earthbox_object::Any,
    obj_dict::Dict{String, T}
)::Nothing where T
    if check_for_earthbox_object(earthbox_object)
        check_and_add!(obj_dict, earthbox_object.name, earthbox_object)
    end
    return nothing
end

function check_for_earthbox_object(
    earthbox_object::Any
)::Bool
    if (
        earthbox_object isa AbstractParameter || 
        earthbox_object isa AbstractEarthBoxArray1D || 
        earthbox_object isa AbstractEarthBoxArray2D
    ) && hasproperty(earthbox_object, :name)
        return true
    else
        return false
    end
end

"""
This function is used to add parameter and array objects to a dictionary. 
Objects are only added if the object name is not already in the dictionary.

# Arguments
- obj_dict::ObjDictType: A dictionary of objects.
- obj_name::String: The name of the object to be added.
- obj::Any: The object parameter or array object to be added to the dictionary.

"""
function check_and_add!(
    obj_dict::ObjDictType,
    obj_name::String, 
    obj::Union{AbstractParameter, AbstractEarthBoxArray1D, AbstractEarthBoxArray2D}
)::Nothing
    if haskey(obj_dict, obj_name)
        error("The name $obj_name is already in the object dictionary. Exiting program.")
    else
        obj_dict[obj_name] = obj
    end
    return nothing
end

"""
Updates an initialization dictionary by converting values to their required types. 
Checks that all required parameters are present in the input dictionary and applies 
type conversion functions to each value.

# Arguments
- input_dict::Dict{String, Union{Float64, Int}}: Input dictionary with parameter values
- required_params::Dict{String, Function}: Dictionary mapping parameter names to their 
    type conversion functions
- use_lists::Bool = true: If true, converts values to single-element vectors

# Returns
- Dict{String, Any}: Updated dictionary with converted values
"""
function update_initialization_dict(
    input_dict::Dict{String, Union{Float64, Int}},
    required_params::Dict{String, Function},
    use_lists::Bool = true
)::Dict{String, Any}
    updated_dict = Dict{String, Any}()
    if input_dict !== nothing
        required_names = keys(required_params)
        input_names = keys(input_dict)

        for required_name in required_names
            if !(required_name in input_names)
                error("initialization_params dict is missing parameter $required_name")
            end
            type_func = required_params[required_name]
            value = input_dict[required_name]
            if use_lists
                updated_dict[required_name] = [type_func(value)]
            else
                updated_dict[required_name] = type_func(value)
            end
        end
    end

    return updated_dict
end

end # end module
