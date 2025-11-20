module MarkerCoordinates

include("options/Options.jl")
include("init_manager/InitManager.jl")

import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.EarthBoxDtypes: InputDictType
import EarthBox.ModelDataContainer: ModelData
import EarthBox: InitializationTools
import EarthBox: OptionTools
import .Options: get_options
import .Options: option_ids
import .Options: option_names
import .InitManager

export initialize!

const COOR_OPTIONS = get_options()

"""
    initialize!(
        model::ModelData;
        marker_distribution::Union{Int, String, Symbol, Nothing}=nothing
    )::Nothing

Initialize marker coordinates.

# Arguments
- `model::ModelData`: The model data object containing model parameters and arrays.
- `marker_distribution::Union{Int, String, Symbol, Nothing}`: Controls the type of initial 
    marker distribution. The marker distribution type is stored in the model data container 
    as an integer ID (`iuse_random`). If `marker_distribution` is nothing the current 
    marker distribution type defined in the model data container will be used. 
    
Available options for `marker_distribution` are as follows:
- `"Regular"`, `:Regular`, or `0`
    - $(COOR_OPTIONS[option_ids[option_names.Regular]].description)
- `"Randomized"`, `:Randomized`, or `1`
    - $(COOR_OPTIONS[option_ids[option_names.Randomized]].description)

Marker distribution type can be accessed from the model data container as follows:
- `iuse_random = model.markers.parameters.distribution.iuse_random.value`

# Returns
- `Nothing`

"""
function initialize!(
    model::ModelData;
    marker_distribution::Union{Int, String, Symbol, Nothing}=nothing
)::Nothing
    option_id = InitializationTools.update_option_id_using_input_option_name(
        get_options(), marker_distribution, model,
        get_option_id_from_model, update_option_id
        )
    option_name = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    @timeit_memit "Finished initializing marker coordinates" begin
        InitManager.initialize!(model, Val(option_name))
    end
    return nothing
end

function get_option_id_from_model(model::ModelData)::Int
    return model.markers.parameters.distribution.iuse_random.value
end

function update_option_id(model::ModelData, option_id::Int)
    model.markers.parameters.distribution.iuse_random.value = option_id
end

function print_option(model::ModelData)
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Marker Coordinates Option")
end

end # module 