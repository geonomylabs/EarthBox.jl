module MarkerTemperature

include("options/Options.jl")
include("parameters/TempInitParameters.jl")
include("init_manager/InitManager.jl")

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.EarthBoxDtypes: InputDictType
import EarthBox: InitializationTools
import EarthBox: OptionTools
import EarthBox.ConversionFuncs: kelvin_to_celsius
import .TempInitParameters: load_parameters!
import .TempInitParameters: FractureZoneKeyNames
import .TempInitParameters: TemperatureWaveKeyNames
import .TempInitParameters: BoxConvectionKeyNames
import .TempInitParameters: LinearKeyNames
import .TempInitParameters: FourLinearSegmentsKeyNames
import .TempInitParameters: AnalyticalThreeLayerKeyNames
import .Options: option_names
import .Options: option_ids
import .Options: get_options

export initialize!

const TEMP_OPTIONS = get_options()
const PDATA = get_eb_parameters()

function make_initial_temperature_model_types_string()::String
    initial_temperature_model_types_string = ""
    for (option_id, option_state) in TEMP_OPTIONS
        option_name = Symbol(option_state.option_name)
        required_parameters = option_state.required_parameters
        required_geometries = option_state.required_geometries
        list_rparams = join(["    - $param" for param in required_parameters], "\n")
        list_rgeometries = join(["    - $geometry" for geometry in required_geometries], "\n")
        initial_temperature_model_types_string *= """
        ## $(option_state.option_name)
        - `initial_temperature_model` **value**: `"$(option_state.option_name)"`, `:$(option_name)`, or $(option_id)
        - **Description**: $(option_state.description)
        - **Required Parameters** (Set using `parameters` dictionary keyword argument):
        $list_rparams
        - **Required Geometries**:
        $list_rgeometries
        """
    end
    return initial_temperature_model_types_string
end

"""
    initialize!(
        model::ModelData; 
        initial_temperature_model::Union{Int, String, Symbol, Nothing}=nothing, 
        kwargs...
    )::Nothing

Initialize the initial marker temperature model. 

# Arguments
- `model::ModelData`
    - Model data container.
- `initial_temperature_model::Union{Int, String, Symbol, Nothing}=nothing`
    - Controls the type of initial temperature model used to initialize marker temperature. 
       See the **Initial Marker Temperature Models** section below for information on available 
       initial temperature models. Required parameters depend on the selected initial 
       temperature model and are set using the parameters dictionary or by using a 
       yamal formatted model input file. The initial temperature model type is stored 
       in the model data container as an integer ID (`itype_temp`) and a corresponding 
       string name (`stype_temp`). If `initial_temperature_model` is nothing the current 
       initial temperature model defined in the model data container will be used. The 
       initial temperature model parameters can be accessed from the model data container as follows:
        - `itype_temp = model.heat_equation.parameters.initial_condition.itype_temp.value`
        - `stype_temp = model.heat_equation.parameters.initial_condition.stype_temp.value`

# Keyword Arguments
- `temperature_uniform::Union{Float64, Nothing}`
    - $(PDATA.temperature_uniform.description)

- `adiabatic_gradient::Union{Float64, Nothing}`
    - $(PDATA.adiabatic_gradient.description)

- `parameters::Union{InputDictType, Nothing}` 
    - Dictionary of required parameters for the selected `initial_temperature_model` model. 
       See the section **Initial Marker Temperature Models** below for information on available 
       initial temperature models and required parameters.
    
# Initial Marker Temperature Models

$(make_initial_temperature_model_types_string())

"""
function initialize!(
    model::ModelData;
    initial_temperature_model::Union{Int, String, Symbol, Nothing}=nothing,
    kwargs...
)::Nothing
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    option_id = InitializationTools.update_option_id_using_input_option_name(
        get_options(), initial_temperature_model, model,
        get_option_id_from_model, update_option_id
        )
    option_symbol = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    load_parameters!(model, option_symbol; kwargs...)
    has_parameters = get(kwargs, :parameters, nothing) !== nothing
    update_parameters_for_transient_cooling!(model, has_parameters)
    InitManager.initialize!(model, Val(option_symbol))
    return nothing
end

""" Set the temperature at the base of the lithosphere using transient cooling parameters.

This is necessary to maintain consistency between initial marker temperature
and transient cooling model. The temperature at the base of the lithosphere
is set to the initial warmer temperature at the base of the lithosphere.
previously calculated during the initialization of the transient cooling model.

"""
function update_parameters_for_transient_cooling!(
    model::ModelData,
    has_parameters::Bool
)::Nothing
    if using_analytical_three_layer(model) && using_bottom_transient(model)
        if has_parameters
            (
                temperature_base_lith_warmer_initial
            ) = model.bcs.parameters.temperature.temperature_base_lith_warmer_initial.value
            model.obj_dict["temperature_base_lith"].value = temperature_base_lith_warmer_initial
            print_info("Found parameters dictionary for analytical three layer model while using transient cooling", level=1)
            print_info("temperature_base_lith_warmer_initial (C): $(kelvin_to_celsius(temperature_base_lith_warmer_initial))", level=2)
            print_info("Set Model Data", level=1)
            print_info("temperature_base_lith (C) (Used by AnalyticalThreeLayer): $(kelvin_to_celsius(model.obj_dict["temperature_base_lith"].value))", level=2)
        end
    end
    return nothing
end

function using_analytical_three_layer(model::ModelData)::Bool
    option_id = get_option_id_from_model(model)
    option_symbol = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    return option_symbol == option_names.AnalyticalThreeLayer
end

function get_option_id_from_model(model::ModelData)::Int
    return model.heat_equation.parameters.initial_condition.itype_temp.value
end

function get_stype_from_model(model::ModelData)::String
    return model.heat_equation.parameters.initial_condition.stype_temp.value
end

function update_option_id(model::ModelData, option_id::Int)::Nothing
    model.heat_equation.parameters.initial_condition.itype_temp.value = option_id
    return nothing
end

function print_option(model::ModelData)::Nothing
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Marker Temperature Option")
    return nothing
end

function using_bottom_transient(model::ModelData)::Bool
    return model.bcs.parameters.temperature.iuse_bottom_transient.value == 1
end

end # module MarkerTemperature 