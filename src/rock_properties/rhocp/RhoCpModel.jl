module RhoCpModel

include("options/Options.jl")
include("update/UpdateManager.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox: InitializationTools
import EarthBox: OptionTools
import .Options: option_names
import .Options: option_ids
import .Options: get_options
import .UpdateManager

function make_rhocp_models_string()::String
    models_string = ""
    for (option_id, option_state) in get_options()
        option_name = Symbol(option_state.option_name)
        models_string *= """
## $(option_state.option_name)
- `rhocp_model` **value**: `"$(option_state.option_name)"`, `:$(option_name)`, or $(option_id)
- **Description**: $(option_state.description)

"""
    end    
    return models_string
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize marker rhocp model.

# Arguments
- `model::ModelData`: Model data container.

- `rhocp_model::Union{Int, String, Symbol, Nothing}=nothing`: 
    RhoCp option id, string name, or symbol. If nothing, the current 
    option id in the model is used either from default value or values
    provided via input files.

"""
function initialize!(
    model::ModelData;
    rhocp_model::Union{Int, String, Symbol, Nothing}=nothing,
)::Nothing
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    _ = InitializationTools.update_option_id_using_input_option_name(
        get_options(), rhocp_model, model,
        get_option_id_from_model, update_option_id
        )
    return nothing
end

function update_rhocp!(model::ModelData, inside_flags::Vector{Int8})
    option_id = get_option_id_from_model(model)
    option_symbol = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    UpdateManager.update_rhocp!(model, inside_flags, Val(option_symbol))
end

function get_option_id_from_model(model::ModelData)::Int
    return model.heat_equation.parameters.rhocp.itype_rhocp.value
end

function get_stype_from_model(model::ModelData)::String
    return model.heat_equation.parameters.rhocp.stype_rhocp.value
end

function update_option_id(model::ModelData, option_id::Int)::Nothing
    model.heat_equation.parameters.rhocp.itype_rhocp.value = option_id
    return nothing
end

function print_option(model::ModelData)::Nothing
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Marker RhoCp Option")
    return nothing
end

end # module RhoCpModel 