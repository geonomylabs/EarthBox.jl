"""
    StaggeredGrid

Module responsible for initializing the basic and staggered grids used for 
discretizing the Stoke-continuity and heat equations.
"""
module StaggeredGrid

include("options/Options.jl")
include("utils/Parameters.jl")
include("utils/Dimensions.jl")
include("utils/Spacing.jl")
include("utils/PressureGridTools.jl")
include("utils/VxGridTools.jl")
include("utils/VyGridTools.jl")
include("utils/VzGridTools.jl")
include("utils/SxyGridTools.jl")
include("utils/SxzGridTools.jl")
include("utils/SyzGridTools.jl")
include("update/UpdateStaggeredGrids.jl")
include("update/UpdateBasicAndStaggeredGrids.jl")
include("types_manager/TypesManager.jl")
include("init_manager/InitManager.jl")

import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.EarthBoxDtypes: InputDictType
import EarthBox.ModelDataContainer: ModelData
import EarthBox: InitializationTools
import EarthBox: OptionTools
import .Options: OptionState
import .Options: get_options
import .Options: option_ids
import .Options: option_names
import .InitManager
import .Parameters

export initialize!

const GRID_OPTIONS = get_options()

"""
    initialize!(
        model::ModelData; 
        grid_type::Union{Int, String, Symbol, Nothing}=nothing, 
        parameters::Union{Dict{String, <:Union{Float64, Int64, String}}, Nothing}=nothing
    )::Nothing

Initialize staggered grid with the specified grid type and parameters. 

# Arguments
- `model::`[ModelData](@ref EarthBox.ModelDataContainer.ModelData)
    - The model data object containing model parameters and arrays. 
- `grid_type`
    - Controls the type of grid. The `grid_type` argument can be an 
       integer ID, a string, a symbol, or `nothing`. Grid type is stored in the model 
       data container as an integer ID (`itype_grid`) and a corresponding string name 
       (`stype_grid`). If `grid_type` is nothing the current grid type defined in the 
       model data container will be used. Grid type parameters can be accessed
       from the model data container as follows:
        - `itype_grid = model.grids.parameters.grid_options.itype_grid.value`
        - `stype_grid = model.grids.parameters.grid_options.stype_grid.value`
- `parameters`
    - Dictionary of parameters used to define the grid type. Currently
       this is only used for `grid_type = :TtypeRefinedGrid`. See the section below 
       called **Parameter Dictionary Keys** for more information.

# Returns
- `Nothing`

# Grid Types
Available grid types (string, symbol, or integer ID) are as follows:
- `UniformGrid`
    - `grid_type` value: `"UniformGrid"`, `:UniformGrid`, or `0`
    - $(GRID_OPTIONS[option_ids[option_names.UniformGrid]].description)
- `TtypeRefinedGrid`
    - `grid_type` value: `"TtypeRefinedGrid"`, `:TtypeRefinedGrid`, or `1` 
    - $(GRID_OPTIONS[option_ids[option_names.TtypeRefinedGrid]].description)
- `Grid51x51UniformRefinementAlongSides`
    - `grid_type` value: `"Grid51x51UniformRefinementAlongSides"`, `:Grid51x51UniformRefinementAlongSides`, or `2` 
    - $(GRID_OPTIONS[option_ids[option_names.Grid51x51UniformRefinementAlongSides]].description)
- `Grid51x51SmoothedRefinementAlongSides`
    - `grid_type` value: `"Grid51x51SmoothedRefinementAlongSides"`, `:Grid51x51SmoothedRefinementAlongSides`, or `3` 
    - $(GRID_OPTIONS[option_ids[option_names.Grid51x51SmoothedRefinementAlongSides]].description)

# Parameter Dictionary Keys

Required parameter keys for grid type `:TtypeRefinedGrid`:
- `"xo_highres"`: Starting x-coordinate of high resolution region in meters
- `"xf_highres"`: Ending x-coordinate of high resolution region in meters
- `"yf_highres"`: Ending y-coordinate of high resolution region in meters
- `"dx_highres"`: Grid spacing in x-direction for high resolution region in meters
- `"dy_highres"`: Grid spacing in y-direction for high resolution region in meters

Optional parameter keys for grid type `:TtypeRefinedGrid`:
- `"iuse_trench"`: Flag to enable trench-based refinement (0 or 1)
- `"iuse_refinement_delay"`: Flag to enable delayed refinement (0 or 1)
- `"refinement_time"`: Time in Myr when refinement should occur

Notes:
- T-type parameters do not need to be provided as input if they were previously defined when 
   initializing the model using [`EarthBoxState`](@ref EarthBoxState).

"""
function initialize!(
    model::ModelData;
    grid_type::Union{Int, String, Symbol, Nothing}=nothing,
    parameters::Union{Dict{String, <:Union{Float64, Int64, String}}, Nothing}=nothing
)::Nothing
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    option_id = InitializationTools.update_option_id_using_input_option_name(
        get_options(), grid_type, model, get_option_id_from_model, update_option_id)
    option_name = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    Parameters.set_parameters!(model, Symbol(option_name), parameters)
    @timeit_memit "Finished initializing staggered grid" begin
        InitManager.initialize!(model, Val(option_name))
    end
    return nothing
end

function get_active_option(
    model::ModelData
)::OptionState
    option_id = get_option_id_from_model(model)
    return get_options()[option_id]
end

function get_boolean_options(model::ModelData)::Dict{Symbol, Bool}
    option_id = get_option_id_from_model(model)
    bools = get_options()[option_id].bools
    return bools
end

function get_stype_from_model(model::ModelData)::String
    return model.grids.parameters.grid_options.stype_grid.value
end

function get_option_id_from_model(model::ModelData)::Int
    return model.grids.parameters.grid_options.itype_grid.value
end

function update_option_id(model::ModelData, option_id::Int)
    model.grids.parameters.grid_options.itype_grid.value = option_id
end

function get_option_name(model::ModelData)::Symbol
    option_id = get_option_id_from_model(model)
    option_name = OptionTools.get_option_symbol_from_id(option_ids, option_id)
    return option_name
end

function print_option(model::ModelData)
    option_id = get_option_id_from_model(model)
    OptionTools.print_option(get_options()[option_id], option_id, "Grid Option")
end

function check_for_ttype_refinement(model::ModelData)::Bool
    return Symbol(get_option_name(model)) == option_names.TtypeRefinedGrid
end

end # module
