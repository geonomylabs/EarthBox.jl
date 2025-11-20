module MarkerPlasticity

include("options/Options.jl")
include("update/UpdateManager.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox: InitializationTools
import EarthBox: OptionTools
import .Options: get_options
import .Options: option_ids
import .Options: option_names

function initialize!(
    model::ModelData;
    marker_plasticity_model::Union{Int, String, Symbol, Nothing}=nothing,
)::Nothing
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    option_id = InitializationTools.update_option_id_using_input_option_name(
        get_options(), marker_plasticity_model, model,
        get_option_id_from_model, update_option_id
        )
    return nothing
end

function update_plastic_yielding!(
    model::ModelData, 
    inside_flags::Vector{Int8}, 
    no_yielding_in_mobile_wall::Bool
)::Nothing
    option_id = get_option_id_from_model(model)
    option_symbol = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    UpdateManager.update_plastic_yielding!(
        model, inside_flags, no_yielding_in_mobile_wall, Val(option_symbol))
    return nothing
end

function get_active_option(
    model::ModelData
)::OptionState
    option_id = get_option_id_from_model(model)
    return get_options()[option_id]
end

function get_option_id_from_model(model::ModelData)::Int
    return model.materials.parameters.material_description.itype_plasticity.value
end

function get_stype_from_model(model::ModelData)::String
    return model.materials.parameters.material_description.stype_plasticity.value
end

function update_option_id(model::ModelData, option_id::Int)
    model.materials.parameters.material_description.itype_plasticity.value = option_id
end

function print_option(model::ModelData)
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Marker Coordinates Option")
end

end # module
