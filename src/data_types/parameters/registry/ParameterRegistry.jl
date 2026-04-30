"""
    ParameterRegistry

Global parameter registry containing all parameters from all model subsystems.

"""
module ParameterRegistry

import EarthBox.TupleTools: merge_named_tuples
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.ParamMacros: @params
import EarthBox.EarthBoxDtypes: InputParamDictType
import EarthBox.EarthBoxDtypes: InputMaterialParamDictType

# --- Section files ---------------------------------------------------------
# Each file under sections/ defines exactly one get_<subsystem>_parameters()
# function. The order below mirrors the merge order in get_eb_parameters().
include("sections/benchmark.jl")
include("sections/carbonate.jl")
include("sections/conversion.jl")
include("sections/gravity.jl")
include("sections/boundary_condition.jl")
include("sections/geometry.jl")
include("sections/grid.jl")
include("sections/heat_equation.jl")
include("sections/interpolation.jl")
include("sections/markers.jl")
include("sections/materials.jl")
include("sections/melting.jl")
include("sections/stokes_continuity.jl")
include("sections/timestep.jl")
include("sections/topography.jl")
include("sections/material_model.jl")
include("sections/material_override.jl")
include("sections/other.jl")


function get_eb_parameters()::NamedTuple
    return merge_named_tuples(
        get_benchmark_parameters(),
        get_carbonate_parameters(),
        get_conversion_parameters(),
        get_gravity_parameters(),
        get_boundary_condition_parameters(),
        get_geometry_parameters(),
        get_grid_parameters(),
        get_heat_equation_parameters(),
        get_interpolation_parameters(),
        get_markers_parameters(),
        get_materials_parameters(),
        get_melting_parameters(),
        get_stokes_continuity_parameters(),
        get_timestep_parameters(),
        get_topography_parameters(),
        get_material_model_parameters(),
        get_material_override_parameters(),
        get_other_parameters()
    )
end

""" Make a dictionary where keys are parameter names and values are Parameter objects.
"""
function get_parameter_dict()::Dict{String, Union{ParameterFloat, ParameterInt, ParameterStr}}
    param_dict = Dict{String, Union{ParameterFloat, ParameterInt, ParameterStr}}()
    eb_parameters = get_eb_parameters()
    for (field_name, param) in pairs(eb_parameters)
        param_dict[param.name] = param
    end
    return param_dict
end

"""
    get_input_parameter_name_list(),

Return a tuple containing:
1. A vector of parameter object names
2. A dictionary mapping parameter names to their type information

Parameter names from input files are checked against this list.

Each dictionary key is a model parameter. Each value is a tuple containing:
(conversion_function, default_value),

where conversion_function is a Julia function that converts string values to the 
appropriate type and the default value is the value used if input is not provided.

Note that the second element of the tuple is the default value for the
parameter. However, this second element is no longer used.
"""
function get_input_parameter_name_list()::Tuple{Vector{String}, InputParamDictType}
    # Create dictionary by mapping each field to a tuple of (parser, value, units)
    eb_parameters = get_eb_parameters()
    input_param_dict = InputParamDictType(p.name => (p.parser, p.value, p.units) for (_, p) in pairs(eb_parameters))
    return collect(keys(input_param_dict)), input_param_dict
end

"""
    check_against_master(param_name::String),

Check input parameter name against master list.
Returns the conversion function for the parameter if it exists.
Throws an error if the parameter name is not valid.
"""
function check_against_master(param_name::String)::Function
    input_param_names, input_param_dict = get_input_parameter_name_list()
    if !(param_name in input_param_names)
        throw(ArgumentError("$param_name is not a valid input parameter name"))
    end
    conv_func = input_param_dict[param_name][1]
    return conv_func
end

""" Get list of missing parameter names that are in master list but not in param_dict and not obsolete.
"""
function get_missing_names(
    param_dict::InputParamDictType,
    master_param_names::Vector{String},
    obsolete_param_names::Vector{String}
)
    missing_names = String[]
    for param_name in master_param_names
        if !(param_name in keys(param_dict)) && !(param_name in obsolete_param_names)
            push!(missing_names, param_name)
        end
    end
    return missing_names
end

function get_master_material_parameters()::Tuple{Vector{String}, InputMaterialParamDictType}
    eb_parameters = get_material_model_parameters()
    master_parameters = InputMaterialParamDictType(
        p.name => (p.parser, p.value, p.units) for (_, p) in pairs(eb_parameters)
    )
    master_names = collect(keys(master_parameters))
    return master_names, master_parameters
end

"""
    check_against_material_master(param_name::String)

Check input parameter name against material master list.
Returns the conversion function for the parameter if it exists.
Throws an error if the parameter name is not valid.
"""
function check_against_material_master(param_name::String)::Function
    master_names, master_parameters = get_master_material_parameters()
    if !(param_name in master_names)
        throw(ArgumentError("$param_name is not a valid material parameter name"))
    end
    conv_func = master_parameters[param_name][1]
    return conv_func
end

""" Get list of missing parameter names that are in master list but not in param_dict and not obsolete.
"""
function get_missing_material_names(
    param_dict::InputMaterialParamDictType,
    master_param_names::Vector{String},
    obsolete_param_names::Vector{String}
)
    missing_names = String[]
    for param_name in master_param_names
        if !(param_name in keys(param_dict)) && !(param_name in obsolete_param_names)
            push!(missing_names, param_name)
        end
    end
    return missing_names
end


end # module