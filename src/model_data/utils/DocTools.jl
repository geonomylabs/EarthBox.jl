module DocTools

import EarthBox.EarthBoxDtypes: AbstractParameterGroup
import EarthBox.EarthBoxDtypes: AbstractArrayGroup
import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import EarthBox.EarthBoxDtypes: AbstractParameterCollection

function make_collection_doc(;
    collection_name::String,
    parameters::Union{AbstractParameterCollection, Nothing} = nothing,
    arrays::Union{AbstractArrayCollection, Nothing} = nothing
)::String
    if parameters isa AbstractParameterCollection && arrays isa Nothing
        return (
            """
            # Accessing Parameters
            $(make_parameter_collection_doc(collection_name=collection_name, parameters=parameters))

            """
        )
    end
    if arrays isa AbstractArrayCollection && parameters isa Nothing
        return (
            """
            # Accessing Arrays
            $(make_array_collection_doc(collection_name=collection_name, arrays=arrays))

            """
        )
    end
    if parameters isa AbstractParameterCollection && arrays isa AbstractArrayCollection
        return (
            """
            # Accessing Parameters
            $(make_parameter_collection_doc(collection_name=collection_name, parameters=parameters))

            # Accessing Arrays
            $(make_array_collection_doc(collection_name=collection_name, arrays=arrays))

            """
        )
    end
end

function make_array_collection_doc(;
    collection_name::String,
    arrays::Union{AbstractArrayCollection, Nothing} = nothing
)::String
    if arrays isa AbstractArrayCollection
        return (
            """
            $(join(
                [get_array_group_doc(collection_name, getfield(arrays, field)) for field in fieldnames(typeof(arrays))], "\n"
            ))
            """
        )
    end
    return ""
end

function get_array_group_doc(
    collection_name::String,
    array_data::AbstractArrayGroup
)::String
    return (
        "#### `model.$(collection_name).arrays.$(array_data.group_name)`\n" *
        join(
            [
                get_array_str(array_data, field)
                for field in fieldnames(typeof(array_data)) if field ∉ [:group_name, :root_dot_path, :title, :obj_list]
            ], "\n"
        )
    )
end

function make_parameter_collection_doc(;
    collection_name::String,
    parameters::Union{AbstractParameterCollection, Nothing} = nothing
)::String
    if parameters isa AbstractParameterCollection
        return (
            """
            $(join(
                [get_parameter_group_doc(collection_name, getfield(parameters, field)) for field in fieldnames(typeof(parameters))], "\n"
            ))
            """
        )
    end
    return ""
end

function get_parameter_group_doc(
    collection_name::String,
    param_data::AbstractParameterGroup
)::String
    return (
        "#### `model.$(collection_name).parameters.$(param_data.group_name)`\n" *
        join(
            [
                get_param_str(param_data, field)
                for field in fieldnames(typeof(param_data)) if field ∉ [:group_name, :root_dot_path, :title, :obj_list]
            ], "\n"
        )
    )
end

function get_param_str(data_obj, field)::String
    return (
        "- `$(getfield(data_obj, field).name).value::$(typeof(getfield(data_obj, field).value))` "
        *": $(getfield(data_obj, field).description) " 
        *": units = $(hasfield(typeof(getfield(data_obj, field)), :units) ? getfield(data_obj, field).units : "None")"
    )
end

function get_array_str(data_obj, field)::String
    return (
        "- `$(getfield(data_obj, field).name).array::$(typeof(getfield(data_obj, field).array))` "
        *": $(getfield(data_obj, field).description) " 
        *": units = $(hasfield(typeof(getfield(data_obj, field)), :units) ? getfield(data_obj, field).units : "None")"
    )
end


end # module 