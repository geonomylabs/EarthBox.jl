module BoundaryConditions

include("bc_calculators/TransientBottomCalculator.jl")
include("options/Options.jl")
include("utils/TransferBC.jl")
include("utils/UpdateBottomVelocityBC.jl")
include("boundary_velocity/BoundaryVelocity.jl")
include("set_boundary_conditions/SetBoundaryConditions.jl")
include("bc_reset/BCResetManager.jl")
include("marker_recycle/MarkerRecycle.jl")
include("bc_managers/Pressure.jl")
include("bc_managers/VelocityFromStrainRate.jl")
include("bc_managers/Temperature.jl")
include("bc_managers/Velocity.jl")
include("bc_managers/VelocityStop.jl")
include("bc_managers/VelocityStep.jl")
include("bc_managers/TransientBottomTemperature.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox: InitializationTools
import EarthBox: OptionTools
import .TransferBC
import .Options: get_options
import .Options: option_ids
import .Options: option_names
import .VelocityStop
import .VelocityStep
import .TransientBottomTemperature

export initialize!, Pressure, VelocityFromStrainRate, Temperature, Velocity, VelocityStop, 
    VelocityStep, TransientBottomTemperature

const BC_OPTIONS = get_options()

function make_bc_model_types_string()::String
    bc_model_types_string = ""
    
    for (option_id, option_state) in BC_OPTIONS
        option_name = Symbol(option_state.option_name)
        bc_model_types_string *= """
## $(option_state.option_name)
- `model_type` **value**: `"$(option_state.option_name)"`, `:$(option_name)`, or $(option_id)
- **Description**: $(option_state.description)

"""
    end
    
    return bc_model_types_string
end

"""
    initialize!(
        model::ModelData;
        model_type::Union{Int, String, Symbol, Nothing}=nothing
    )::Nothing

Initialize integrated boundary condition model type for coupled Stokes-continuity
and heat equations.

# Arguments
- `model::`[`ModelData`](@ref ModelData): The model data object containing model parameters and arrays.
- `model_type::Union{Int, String, Symbol, Nothing}`: Controls the type of integrated boundary 
    condition model type for the Stokes-continuity and heat equations. See the **Integrated Boundary 
    Condition Model Types** section below for information on available boundary condition model types. 
    The boundary condition model type is stored in the model data container as an integer ID (`itype_bc`) 
    and a corresponding string name (`stype_bc`). If `model_type` is nothing the current boundary condition 
    type defined in the model data container will be used. The boundary condition model type parameters can be 
    accessed from the model data container as follows:
    - `itype_bc = model.bcs.parameters.bc_options.itype_bc.value`
    - `stype_bc = model.bcs.parameters.bc_options.stype_bc.value`

# Returns
- `Nothing`

# Integrated Boundary Condition Model Types

$(make_bc_model_types_string())

"""
function initialize!(
    model::ModelData;
    model_type::Union{Int, String, Symbol, Nothing}=nothing,
)::Nothing
    # Ensure that the option id (i.e. itype_bc) is consistent with the
    # input boundary condition model type (i.e. stype_bc)
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    _ = InitializationTools.update_option_id_using_input_option_name(
        get_options(), model_type, model,
        get_option_id_from_model, update_option_id
    )
    return nothing
end

function print_option(model::ModelData)::Nothing
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Boundary Conditions Option")
    return nothing
end

"""
    update_transient_boundary_conditions!(model::ModelData)::Nothing

Update transient boundary conditions at the start of each time step for velocity stepping,
velocity stopping and transient bottom temperature.

"""
function update_transient_boundary_conditions!(model::ModelData)::Nothing
    @timeit_memit "Finished updating transient boundary conditions" begin
        reset_for_velocity_step!(model)
        reset_for_velocity_stop!(model)
        reset_bc!(model)
        TransientBottomTemperature.update_bottom_temperature!(model)
    end
    return nothing
end

"""
    reset_for_velocity_step!(model::ModelData)::Nothing

Reset boundary conditions for velocity stepping. If a velocity step is used, recalculate 
conservative upper and lower velocity boundary conditions and reset boundary condition
coefficient arrays for Stokes equations and heat equation.
"""
function reset_for_velocity_step!(model::ModelData)::Nothing
    VelocityStep.reset_bc_for_veloc_step!(model)
    return nothing
end

"""
    reset_for_velocity_stop!(model::ModelData)::Nothing

Set velocity boundary conditions to zero.
"""
function reset_for_velocity_stop!(model::ModelData)::Nothing
    VelocityStop.reset_bc_for_veloc_stop!(model)
    return nothing
end

function reset_bc!(model::ModelData)::Nothing
    option_id = get_option_id_from_model(model)
    option_name = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    BCResetManager.reset_bcs!(model, Val(option_name))
    TransferBC.transfer_velocity_bcs_to_component_arrays!(model)
    return nothing
end

function manage_marker_recycling!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "Finished marker recycling" begin    
        option_id = get_option_id_from_model(model)
        option_name = OptionTools.get_option_symbol_from_id(option_ids, option_id)
        MarkerRecycle.recycle_markers!(model, inside_flags, Val(option_name))
    end
    return nothing
end

function get_boolean_options(model::ModelData)::Dict{Symbol, Bool}
    option_id = get_option_id_from_model(model)
    bools = get_options()[option_id].bools
    return bools
end

function get_option_id_from_model(model::ModelData)::Int
    return model.bcs.parameters.bc_options.itype_bc.value
end

function get_stype_from_model(model::ModelData)::String
    return model.bcs.parameters.bc_options.stype_bc.value
end

function update_option_id(model::ModelData, option_id::Int)::Nothing
    model.bcs.parameters.bc_options.itype_bc.value = option_id
    return nothing
end

function get_option_name(model::ModelData)::Symbol
    option_id = get_option_id_from_model(model)
    option_name = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    return option_name
end

end # module 