module OptionTools

import DataStructures: OrderedDict
import EarthBox: PrintFuncs
import EarthBox.EarthBoxDtypes: AbstractOptionType

export print_option, print_all, check_option_id, get_id, get_name, 
    get_option_symbol_from_id

function make_option_names(option_ids::Dict{Symbol, Int})::NamedTuple
    return NamedTuple(keys(option_ids) .=> keys(option_ids))
end

struct OptionState <: AbstractOptionType
    option_name::Union{String, Nothing}
    description::Union{String, Nothing}
    bools::Union{Dict{Symbol, Bool}, Nothing}
    required_domains::Union{Vector{String}, Nothing}
    required_types::Union{Vector{String}, Nothing}
    required_parameters::Union{Vector{String}, Nothing}
    required_geometries::Union{Vector{String}, Nothing}
    time_steps::Union{Vector{Int}, Nothing} # Time steps to post-process (for benchmarks)
    y_index::Union{Int, Nothing} # Y-index to post-process (for benchmarks)

    function OptionState(;
        option_name::Union{String, Nothing} = nothing,
        description::Union{String, Nothing} = nothing,
        bools::Union{Dict{Symbol, Bool}, Nothing} = nothing,
        required_domains::Union{Vector{String}, Nothing} = nothing,
        required_types::Union{Vector{String}, Nothing} = nothing,
        required_parameters::Union{Vector{String}, Nothing} = nothing,
        required_geometries::Union{Vector{String}, Nothing} = nothing,
        time_steps::Union{Vector{Int}, Nothing} = nothing,
        y_index::Union{Int, Nothing} = nothing,
    )
        if bools === nothing
            bools = Dict{Symbol, Bool}()
        end
        new(
            option_name, description, bools, required_domains, 
            required_types, required_parameters, required_geometries, 
            time_steps, y_index
            )
    end
end

function Base.print(opt::OptionState)
    PrintFuncs.print_option_info(
        opt.option_name === nothing ? "None" : opt.option_name,
        opt.description === nothing ? "None" : opt.description
    )
end

function print_option(option::OptionState, option_id::Int, name::String)
    PrintFuncs.print_option_category_info(name)
    PrintFuncs.print_option_id(option_id)
    PrintFuncs.print_option_info(option.option_name, option.description)
end

function print_all(options::Union{Dict{Int, OptionState}, OrderedDict{Int, OptionState}}, name::String)
    for option_id in keys(options)
        print_option(options[option_id], option_id, name)
    end
end

function check_option_id(options::Union{Dict{Int, OptionState}, OrderedDict{Int, OptionState}}, option_id::Int)
    if !haskey(options, option_id)
        throw(ArgumentError("option_id $option_id is not a valid option."))
    end
end

function get_id(options::Union{Dict{Int, OptionState}, OrderedDict{Int, OptionState}}, option_name::Union{String, Symbol})::Int
    if option_name isa Symbol
        option_name = String(option_name)
    end
    for (option_id, opt) in options
        if opt.option_name == option_name
            return option_id
        end
    end
    throw(ArgumentError("Could not find option with name $option_name in options dictionary."))
end

function get_name(options::Union{Dict{Int, OptionState}, OrderedDict{Int, OptionState}}, option_id::Int)::String
    if !haskey(options, option_id)
        throw(ArgumentError("Could not find option with id $option_id in options dictionary."))
    end
    return options[option_id].option_name
end

function get_option_symbol_from_id(
    option_ids::Dict{Symbol, Int},
    option_id::Int
)::Symbol
    for (option_name, id) in option_ids
        if id == option_id
            return option_name
        end
    end
    error("Option id $option_id not found in option_ids")
end

end # module 