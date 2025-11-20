module Sealevel

include("options/Options.jl")
include("relative_base_level/RelativeBaseLevel.jl")
include("update/UpdateSealevel.jl")

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox: InitializationTools
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox: OptionTools
import .Options: get_options
import .Options: option_ids
import .Options: option_names

const PDATA = get_eb_parameters()
const SEALEVEL_OPTIONS = get_options()

struct ValidInputNames
    y_sealevel::Symbol # meters
    base_level_shift::Symbol # meters
    base_level_shift_end_time::Symbol # millions of years
end

function make_sealevel_option_string()::String
    sealevel_option_string = ""
    for (option_id, option_state) in SEALEVEL_OPTIONS
        option_name = Symbol(option_state.option_name)
        sealevel_option_string *= """
## $(option_state.option_name)
- `option_name` **value**: `"$(option_state.option_name)"`, `:$(option_name)`, or $(option_id)
- **Description**: $(option_state.description)
"""
    end
    return sealevel_option_string
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize sealevel model parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData)
    - The model data container containing the model parameters and arrays.
- `option_name::Union{Int, String, Symbol, Nothing}=nothing`
    - Controls the type of sealevel option. See the **Sealevel Options** section below for 
        information on available sealevel options. The sealevel option is stored in the model data 
        container as an integer ID (`itype_sealevel`) and a corresponding string name (`stype_sealevel`). 
        If `option_name` is nothing the current sealevel option defined in the model data container 
        will be used. The sealevel option parameters can be accessed from the model data container as follows:
        - `itype_sealevel = model.topography.parameters.sealevel.itype_sealevel.value`
        - `stype_sealevel = model.topography.parameters.sealevel.stype_sealevel.value`

# Keyword Arguments
- `$(PDATA.y_sealevel.name)::Float64`
    - $(PDATA.y_sealevel.description)
- `$(PDATA.base_level_shift.name)::Float64`
    - $(PDATA.base_level_shift.description)
- `$(PDATA.base_level_shift_end_time.name)::Float64`
    - $(PDATA.base_level_shift_end_time.description)

---
# Sealevel Options
---
$(make_sealevel_option_string())

"""
function initialize!(
    model::ModelData;
    option_name::Union{Int, String, Symbol, Nothing}=nothing,
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    y_sealevel = get(kwargs, :y_sealevel, nothing)
    set_initial_water_level!(model, y_sealevel)
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    _ = InitializationTools.update_option_id_using_input_option_name(
        get_options(), option_name, model, get_option_id_from_model, update_option_id)
    return nothing
end

function set_initial_water_level!(
    model::ModelData,
    y_sealevel::Union{Float64, Nothing}=nothing
)::Nothing
    if y_sealevel !== nothing
        model.topography.parameters.sealevel.y_water_ini.value = y_sealevel
    end
    # ensure that y_sealevel is initialized to y_water_ini
    model.topography.parameters.sealevel.y_sealevel.value = model.topography.parameters.sealevel.y_water_ini.value
    return nothing
end

function get_active_option_symbol(model::ModelData)::Symbol
    option_id = get_option_id_from_model(model)
    option_symbol = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    return option_symbol
end

""" Check if sealevel is defined by the upper left edge of the crust.
"""
function use_upper_left_edge_of_crust_to_define_sealevel(model::ModelData)::Bool
    option_symbol = get_active_option_symbol(model)
    return option_symbol == option_names.LeftEdge
end

function update_option_id(model::ModelData, option_id::Int)::Nothing
    model.topography.parameters.sealevel.itype_sealevel.value = option_id
    return nothing
end

function get_option_id_from_model(model::ModelData)::Int
    return model.topography.parameters.sealevel.itype_sealevel.value
end

function get_stype_from_model(model::ModelData)::String
    return model.topography.parameters.sealevel.stype_sealevel.value
end

function print_option(model::ModelData)
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Sealevel Option")
end

""" Update sea level.

Options include left edge, average pressure and fixed sea level. The
left edge option uses the upper-left corner of the lithosphere,
average pressure uses a reference lithosphere, the average pressure of 
the lithosphere and an assumption of isostasy, and fixed sea level
uses a fixed value for sea level.
"""
function update_sealevel!(
    model::ModelData
)::Nothing
    iuse_topo = model.topography.parameters.model_options.iuse_topo.value
    if iuse_topo == 1
        UpdateSealevel.update_sealevel!(model, Val(get_active_option_symbol(model)))
    end
    return nothing
end

end # module 