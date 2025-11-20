module MarkerFriction

include("utils/FrictionRandomizer.jl")
include("options/Options.jl")
include("init_manager/InitManager.jl")

import Random: rand
import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.EarthBoxDtypes: InputDictType
import EarthBox: InitializationTools
import EarthBox.ParameterGroupTools: set_group_parameters!
import EarthBox.Rheology.PlasticFailure.RandomFrictionAngles:
    randomize_initial_friction_angles!
import EarthBox: OptionTools
import .Options: get_options
import .Options: option_ids
import .Options: option_names
import .InitManager

struct ValidInputNames
    delta_fric_coef::Symbol
    iuse_random_fric_time::Symbol
    randomization_factor::Symbol
end

const FRIC_OPTIONS = get_options()
const PDATA = get_eb_parameters()

function make_friction_initialization_models_string()::String
    friction_initialization_model_string = ""
    for (option_id, option_state) in FRIC_OPTIONS
        option_name = Symbol(option_state.option_name)
        friction_initialization_model_string *= """
        ## $(option_state.option_name)
        - `initialization_model` **value**: `"$(option_state.option_name)"`, `:$(option_name)`, or $(option_id)
        - **Description**: $(option_state.description)
        """
    end
    return friction_initialization_model_string
end

"""
    initialize!(
        model::ModelData; 
        initialization_model::Union{Int, String, Symbol, Nothing}=nothing, kwargs...
    )::Nothing

Manage marker friction coefficients initialization and time-dependent models.

# Arguments
- `model::ModelData`: The model data container containing the model parameters and arrays.
- `initialization_model::Union{Int, String, Symbol, Nothing}`: 
    - Controls the type of marker friction coefficients initialization model.
       initialization model for marker friction coefficients. See the 
       **Marker Friction Coefficient Initialization Models** section below for information on 
       available initialization models. If the `initialization_model` is not provided, the default 
       initialization model will be used. The initialization model is stored in the model data 
       container as an integer flag (`iuse_random_fric`). The initialization model parameter flag 
       can be accessed from the model data container as follows:
        - `iuse_random_fric = model.materials.parameters.random_friction.iuse_random_fric.value`

# Keyword Arguments
- `delta_fric_coef::Union{Float64, Nothing}=nothing`:
    - $(PDATA.delta_fric_coef.description)
- `iuse_random_fric_time::Union{Bool, Nothing}=nothing`:
    - $(PDATA.iuse_random_fric_time.description)
- `randomization_factor::Union{Float64, Nothing}=nothing`:
    - $(PDATA.randomization_factor.description)

# Marker Friction Initialization Models

$(make_friction_initialization_models_string())

"""
function initialize!(
    model::ModelData;
    initialization_model::Union{Int, String, Symbol, Nothing}=nothing,
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    
    option_id = InitializationTools.update_option_id_using_input_option_name(
        get_options(), initialization_model, model,
        get_option_id_from_model, update_option_id
    )
    
    option_name = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    InitManager.initialize!(model, Val(option_name))
    
    return nothing
end

"""
    randomize_marker_initial_friction_angle!(model::ModelData)::Nothing

Manage randomization of initial friction angles for markers.

This method manages the randomization of the initial friction angle 
*θ°ₘ* for each marker. The sine of this randomized angle (i.e. friction
coefficient) is stored in the array `marker_fric_ini`.

# Arguments
- `model::ModelData`: The model data.

# Returns
- `Nothing`: This function does not return anything.
"""
function randomize_marker_initial_friction_angle!(
    model::ModelData, 
)::Nothing
    @timeit_memit "Finished randomizing marker initial friction angle" begin
        option_id = get_option_id_from_model(model)
        if option_id == option_ids[:Randomized]
            randomize_initial_friction_angles!(model)
        end
    end
    return nothing
end

function print_option(model::ModelData)
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Marker Friction Option")
end

function get_option_id_from_model(model::ModelData)::Int
    return model.materials.parameters.random_friction.iuse_random_fric.value
end

function update_option_id(model::ModelData, option_id::Int)::Nothing
    model.materials.parameters.random_friction.iuse_random_fric.value = option_id
    return nothing
end

end # module 