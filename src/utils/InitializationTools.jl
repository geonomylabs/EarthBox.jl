module InitializationTools

import DataStructures: OrderedDict
import EarthBox.ModelDataContainer: ModelData
import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: get_id
import EarthBox.OptionTools: check_option_id

export update_use_flags, update_option_id_using_input_option_name, 
    sync_option_id_with_stype, get_option_id_from_option_name

""" Update use flags.

If buse (boolean use flag) is not Nothing then use it to define iuse
(integer use flag). If buse is Nothing then use iuse to define buse.
"""
function update_use_flags(
    buse::Union{Bool, Nothing},
    iuse::Int
)::Tuple{Bool, Int}
    if buse isa Bool
        iuse = buse ? 1 : 0
    else
        buse = iuse == 1
    end
    return buse, iuse
end

""" 
Translate `option_name` to `option_id` and update `option_id` in model data 
container. If `option_name` is Nothing then the current (initial) `option_id` in 
the model data container `model` is used and returned. Otherwise, if 
`option_name` is a non-Nothing and a new `option_id` is encountered then the 
`option_id` in the model data object is updated and returned.

# Arguments
- `options::Union{Dict{Int, OptionState}, OrderedDict{Int, OptionState}}`
    - Dictionary of option ids and option states.
- `option_name::Union{Int, String, Symbol, Nothing}`
    - Option name to convert to option id.
- `model::ModelData`
    - Model data container.
- `get_option_id_from_model_func::Function`
    - Function to get the option id from the model data container.
- `update_option_id_func::Function`
    - Function to update the option id in the model data container.

# Returns
- `option_id::Int`
    - Option id.
"""
function update_option_id_using_input_option_name(
    options::Union{Dict{Int, OptionState}, OrderedDict{Int, OptionState}},
    option_name::Union{Int, String, Symbol, Nothing},
    model::ModelData,
    get_option_id_from_model_func::Function,
    update_option_id_func::Function
)::Int
    option_id_init = get_option_id_from_model_func(model)
    option_id, update = get_option_id_from_option_name(
        options, option_name, option_id_init)
    if update
        update_option_id_func(model, option_id)
    end
    return option_id
end

""" Sync option id with stype parameter.

Sync option_id (itype parameter) with stype parameter. This function 
essentially forces consistency between the option id and the stype parameter
defined in the model data container. If stype parameter is not defined then 
option id (itype) is not modified.

This method is called with initialization of all managers that define the
main components of the model and allows string names to be used, which are
easier to understand, in input files instead of integer option ids.
"""
function sync_option_id_with_stype(
    options::Union{Dict{Int, OptionState}, OrderedDict{Int, OptionState}},
    model::ModelData,
    get_stype_from_model_func::Function,
    update_option_id_func::Function
)::Nothing
    stype = get_stype_from_model_func(model)
    if stype != "None"
        option_id = get_id(options, stype)
        update_option_id_func(model, option_id)
    end
    return nothing
end

function get_option_id_from_option_name(
    options::Union{Dict{Int, OptionState}, OrderedDict{Int, OptionState}},
    option_name::Union{Int, String, Symbol, Nothing},
    option_id_init::Int
)::Tuple{Int, Bool}
    update = false
    if option_name !== nothing
        option_id = option_name
        if option_name isa Int
            check_option_id(options, option_name)
        elseif option_name isa String || option_name isa Symbol
            option_id = get_id(options, option_name)
        end
        if option_id != option_id_init
            update = true
        end
    else
        option_id = option_id_init
    end
    return option_id, update
end

end # module 