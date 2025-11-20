module GlobalPlasticityLoop

include("options/Options.jl")
include("utils/Initialize.jl")
include("utils/LoopInformation.jl")
include("utils/LoopTermination.jl")
include("utils/UpdateLoopParameters.jl")
include("utils/Convergence.jl")
include("plasticity_loop/PlasticityLoop.jl")

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit, print_info
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.StokesContinuitySolver: interpolate_staggered_grid_velocity_to_basic_grid!
import EarthBox: InitializationTools
import EarthBox.LoopInputStruct: LoopInput
import EarthBox: OptionTools
import .Options: OptionState
import .Options: get_options
import .Options: option_ids
import .Options: option_names
import .PlasticityLoop

const PDATA = get_eb_parameters()
const PLASTIC_LOOP_OPTIONS = get_options()

struct ValidInputNames
    tolerance_picard::Symbol
    nglobal::Symbol
end

function make_global_plasticity_loop_string()::String
    global_plasticity_loop_string = ""
    for (option_id, option_state) in PLASTIC_LOOP_OPTIONS
        option_name = Symbol(option_state.option_name)
        global_plasticity_loop_string *= """
## $(option_state.option_name)
- `global_plasticity_loop` **value**: `"$(option_state.option_name)"`, `:$(option_name)`, or $(option_id)
- **Description**: $(option_state.description)
"""
    end
    return global_plasticity_loop_string
end

"""
    initialize!(
        model::ModelData;
        global_plasticity_loop::Union{Int, Symbol, String, Nothing}=nothing,
        kwargs...
    )::Nothing

Initialize the global plasticity loop (Picard iteration loop]).

# Arguments
- `model::`[`ModelData`](@ref ModelData): The model data container containing 
    model parameters and arrays.

# Keyword Arguments
- `global_plasticity_loop::Union{Int, String, Nothing}=nothing`:
    - Controls the type of global plasticity loop. See the **Global Plasticity Loop Types** section 
       below for information on available global plasticity loop types. If `global_plasticity_loop` 
       is nothing the option id from the model data container will be used to define global plasticity 
       loop type. The global plasticity loop type is stored in the model data container as an integer 
       ID (`itype_global`) and a corresponding string name (`stype_global`). The global plasticity loop 
       type parameters can be accessed from the model data container as follows:
        - `itype_global = model.stokes_continuity.parameters.picard.itype_global.value`
        - `stype_global = model.stokes_continuity.parameters.picard.stype_global.value`
- `$(PDATA.tolerance_picard.name)::Float64`:
    - $(PDATA.tolerance_picard.description)
- `$(PDATA.nglobal.name)::Int64`:
    - $(PDATA.nglobal.description)
---
# Global Plasticity Loop Types
---
$(make_global_plasticity_loop_string())

"""
function initialize!(
    model::ModelData;
    global_plasticity_loop::Union{Int, Symbol, String, Nothing}=nothing,
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    _ = InitializationTools.update_option_id_using_input_option_name(
        get_options(), global_plasticity_loop, model, get_option_id_from_model, update_option_id)
    return nothing
end

function print_option(model::ModelData)
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Global Plasticity Loop Option")
end

"""
Check global plasticity loop option. If option is set to marker_plasticity_loop as opposed to
nodal_plasticity_loop return true.
"""
function is_global_marker_plasticity_loop(model::ModelData)::Bool
    option_id = get_option_id_from_model(model)
    option_symbol = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    return option_symbol == option_names.MarkerPlasticityLoop
end

"""
Solve Stokes-continuity equations using Picard iteration method. This method calls the active 
plasticity loop function to solve the Stokes-continuity equations for velocity and pressure using 
the Picard iteration method. The plasticity loop function is either marker based or node based 
depending on the option selected. The new staggered grid velocity is interpolated to the basic 
grid. Deviatoric grid stress and viscoplastic viscosity on the staggered grid are also updated.
"""
function run_picard_loop!(
    model::ModelData, 
    loop_input::LoopInput,
    inside_flags::Vector{Int8}
)::Nothing
    option_id = get_option_id_from_model(model)
    print_info("Starting picard loop", level=2)
    option_symbol = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    @timeit_memit "# Finished picard loop" begin
        PlasticityLoop.run_picard_loop(model, loop_input, inside_flags, Val(option_symbol))
    end
    @timeit_memit "Finished interpolating staggered grid velocity to basic grid" begin
        interpolate_staggered_grid_velocity_to_basic_grid!(model)
    end
    return nothing
end

function get_option_id_from_model(model::ModelData)::Int
    return model.stokes_continuity.parameters.picard.itype_global.value
end

function get_stype_from_model(model::ModelData)::String
    return model.stokes_continuity.parameters.picard.stype_global.value
end

function update_option_id(model::ModelData, option_id::Int)::Nothing
    model.stokes_continuity.parameters.picard.itype_global.value = option_id
    return nothing
end

end # module 