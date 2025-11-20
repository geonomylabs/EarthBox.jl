module Parameters

import EarthBox.EarthBoxDtypes: AbstractParameter


"""
    ParameterFloat

A mutable struct representing a parameter with a Float64 value.

# Fields
- `value::Float64`: The numerical value of the parameter
- `name::String`: The name of the parameter 
- `units::String`: The units of the parameter
- `description::String`: A description of what the parameter represents
- `parser::Function`: A function to parse the parameter value from a string
"""
mutable struct ParameterFloat <: AbstractParameter
    value::Float64
    name::String
    units::String
    description::String
    parser::Function
    function ParameterFloat(
        value::Union{Float64, Int},
        name::String, 
        units::String, 
        description::String, 
        parser::Function = x -> parse(Float64, x)
    )
        if typeof(value) == Int
            value = convert(Float64, value)
        end
        return new(value, name, units, description, parser)
    end
end

function set_parameter_float!(parameter::ParameterFloat, value::Float64)
    parameter.value = value
end

"""
    ParameterInt

A mutable struct representing a parameter with an Int64 value.

# Fields
- `value::Int64`: The numerical value of the parameter
- `name::String`: The name of the parameter 
- `units::String`: The units of the parameter
- `description::String`: A description of what the parameter represents
- `parser::Function`: A function to parse the parameter value from a string
"""
mutable struct ParameterInt <: AbstractParameter
    value::Int64
    name::String
    units::String
    description::String
    parser::Function
    function ParameterInt(
        value::Union{Int64, Float64},
        name::String, 
        units::String, 
        description::String, 
        parser::Function = x -> parse(Int64, x)
    )
        if typeof(value) == Float64
            value = convert(Int64, floor(value))
        end
        return new(value, name, units, description, parser)
    end
end

function set_parameter_int!(parameter::ParameterInt, value::Int64)
    parameter.value = value
end

"""
    ParameterStr

A mutable struct representing a parameter with a String value.

# Fields
- `value::String`: The string value of the parameter
- `name::String`: The name of the parameter 
- `units::String`: The units of the parameter
- `description::String`: A description of what the parameter represents
- `parser::Function`: A function to parse the parameter value from a string
"""
mutable struct ParameterStr <: AbstractParameter
    value::String
    name::String
    units::String
    description::String
    parser::Function
    function ParameterStr(
        value::String, 
        name::String, 
        units::String, 
        description::String, 
        parser::Function = string
    )
        return new(value, name, units, description, parser)
    end
end

function set_parameter_str!(parameter::ParameterStr, value::String)
    parameter.value = value
end

function print_parameter(
    parameter::Union{ParameterFloat, ParameterInt, ParameterStr}
)::Nothing
    println("\n parameter name: ", parameter.name)
    println(repeat("-", 79))
    println("description: ", parameter.description)
    println("units: ", parameter.units)
    println("value: ", parameter.value)
    return nothing
end

function print_param_obj(param_obj)
    println("$(param_obj.name) : $(param_obj.value) : $(param_obj.units) : $(param_obj.description)")
end

function convert_to_float64(value::Union{Float64, Int64, String, Bool})::Float64
    if typeof(value) == Float64
        return value
    elseif typeof(value) == Int64
        return convert(Float64, value)
    elseif typeof(value) == String
        return parse(Float64, value)
    elseif typeof(value) == Bool
        return value ? 1.0 : 0.0
    end
end

function convert_to_int64(value::Union{Float64, Int64, String, Bool})::Int64
    if typeof(value) == Int64
        return value
    elseif typeof(value) == Float64
        return convert(Int64, value)
    elseif typeof(value) == String
        return floor(Int64, parse(Float64, value))
    elseif typeof(value) == Bool
        return value ? 1 : 0
    end
end

end # module
