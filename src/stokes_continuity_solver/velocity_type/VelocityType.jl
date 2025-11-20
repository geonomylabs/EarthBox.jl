module VelocityType

include("options/Options.jl")

import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ModelDataContainer: ModelData
import EarthBox: InitializationTools
import EarthBox: OptionTools
import .Options: OptionState
import .Options: get_options
import .Options: option_ids
import .Options: option_names

const VEL_TYPE_OPTIONS = get_options()

function make_velocity_type_string()::String
    velocity_type_string = ""
    for (option_id, option_state) in VEL_TYPE_OPTIONS
        option_name = Symbol(option_state.option_name)
        velocity_type_string *= """
## $(option_state.option_name)
- `velocity_type` **value**: `"$(option_state.option_name)"`, `:$(option_name)`, or $(option_id)
- **Description**: $(option_state.description)

"""
    end
    return velocity_type_string
end

"""
    initialize!(
        model::ModelData;
        velocity_type::Union{Int, String, Symbol, Nothing}=nothing
    )::Nothing

Initialize velocity type for the Stokes-continuity solver allowing the user to
turn off the Stokes-continuity solver and prescribe the velocity field using a 
solid body rotation model or set the velocity to zero for testing purposes.

# Arguments
- `model::ModelData`: The model data container containing model parameters and arrays.

# Keyword Arguments
- `velocity_type`: Controls the type of velocity calculation. The `velocity_type` 
  can be an integer ID, a string, a symbol, or nothing. See the **Velocity Types** section below 
  for information on available velocity types. If `velocity_type` is nothing the option id from 
  the model data container will be used to define velocity type. The velocity type is stored in 
  the model data container as an integer ID (`itype_velocity`) and a corresponding string name 
  (`stype_velocity`). The velocity type parameters can be accessed from the model data container as follows:
  - `itype_velocity = model.stokes_continuity.parameters.velocity_calc_option.itype_velocity.value`
  - `stype_velocity = model.stokes_continuity.parameters.velocity_calc_option.stype_velocity.value`

# Returns
- `Nothing`

# Velocity Types
$(make_velocity_type_string())

"""
function initialize!(
    model::ModelData;
    velocity_type::Union{Int, String, Symbol, Nothing}=nothing,
)::Nothing
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    _ = InitializationTools.update_option_id_using_input_option_name(
        get_options(), velocity_type, model, get_option_id_from_model, update_option_id)
    return nothing
end

function is_velocity_from_stokes_solver(model::ModelData)::Bool
    """ Check if velocity is calculated from Stokes-continuity solver.
    """
    option_id = get_option_id_from_model(model)
    option_symbol = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    check = false
    if option_symbol == option_names.VelocityFromStokesSolver
        check = true
    end
    return check
end

function is_velocity_from_solid_body_rotation(model::ModelData)::Bool
    """ Check if velocity is calculated from solid body rotation model.
    """
    option_id = get_option_id_from_model(model)
    option_symbol = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    check = false
    if option_symbol == option_names.KinematicSolidBodyRotation
        check = true
    end
    return check
end

function get_active_option(
    model::ModelData
)::OptionState
    option_id = get_option_id_from_model(model)
    return get_options()[option_id]
end

function get_stype_from_model(model::ModelData)::String
    return model.stokes_continuity.parameters.velocity_calc_option.stype_velocity.value
end

function get_option_id_from_model(model::ModelData)::Int
    return model.stokes_continuity.parameters.velocity_calc_option.itype_velocity.value
end

function update_option_id(model::ModelData, option_id::Int)
    model.stokes_continuity.parameters.velocity_calc_option.itype_velocity.value = option_id
end

function print_option(model::ModelData)
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Velocity Type Option")
end

end # module 