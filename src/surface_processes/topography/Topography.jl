module Topography

include("options/Options.jl")
include("utils/TopographyFuncs.jl")
include("utils/AdvectionTools.jl")
include("utils/CarbonateAdvect.jl")
include("topo_type/TopoTypeManager.jl")

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit, print_info
import EarthBox: InitializationTools
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ModelDataContainer.TopographyContainer.ArrayCollection: initialize_topo_array
import EarthBox.Markers.MarkerMaterials.Options: option_ids as mat_option_ids
import EarthBox.Markers.MarkerMaterials.Options: option_names as mat_option_names
import EarthBox: OptionTools
import .Options: option_ids
import .Options: option_names
import .Options: get_options

export initialize!

const PDATA = get_eb_parameters()
const ADVECTION_OPTIONS = get_options()

function make_node_advection_string()::String
    node_advection_string = ""
    for (option_id, option_state) in ADVECTION_OPTIONS
        option_name = Symbol(option_state.option_name)
        node_advection_string *= """
## $(option_state.option_name)
- `node_advection` **value**: `"$(option_state.option_name)"`, `:$(option_name)`, or $(option_id)
- **Description**: $(option_state.description)
"""
    end
    return node_advection_string
end

struct ValidInputNames
    iuse_topo::Symbol
    nsmooth_top_bottom::Symbol # number of topo nodes used for smoothing
    marker_search_factor::Symbol # search factor for marker search radius
    toponum::Symbol # number of topo nodes
    topo_xsize::Symbol # meters
    dx_topo::Symbol # meters
    erosion_rate::Symbol # uniform erosion rate m/s
    sedimentation_rate::Symbol # uniform sedimentation rate m/s
end


"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize topography model parameters.

# Arguments
- `model::`[`ModelData`](@ref ModelData)
    - The model data container containing the model parameters and arrays.
- `node_advection::Union{Int, String, Symbol, Nothing}`
    - Controls the type of topography node advection scheme. See the 
        **Node Advection Schemes** section below for information on available node advection schemes.
        The node advection scheme is stored in the model data container as an integer ID (`itype_topo`) 
        and a corresponding string name (`stype_topo`). If `node_advection` is nothing the current 
        node advection scheme defined in the model data container will be used. The node advection scheme 
        parameters can be accessed from the model data container as follows:
        - `itype_topo = model.topography.parameters.model_options.itype_topo.value`
        - `stype_topo = model.topography.parameters.model_options.stype_topo.value`

# Keyword Arguments
- `$(PDATA.iuse_topo.name)::Int64`
    - $(PDATA.iuse_topo.description)
- `$(PDATA.nsmooth_top_bottom.name)::Int64`
    - $(PDATA.nsmooth_top_bottom.description)
- `$(PDATA.marker_search_factor.name)::Float64`
    - $(PDATA.marker_search_factor.description)
- `$(PDATA.toponum.name)::Int64`
    - $(PDATA.toponum.description)
- `$(PDATA.topo_xsize.name)::Float64`
    - $(PDATA.topo_xsize.description)
- `$(PDATA.dx_topo.name)::Float64`
    - $(PDATA.dx_topo.description)
- `$(PDATA.erosion_rate.name)::Float64`
    - $(PDATA.erosion_rate.description)
- `$(PDATA.sedimentation_rate.name)::Float64`
    - $(PDATA.sedimentation_rate.description)
---
# Node Advection Schemes
---
$(make_node_advection_string())

"""
function initialize!(
    model::ModelData;
    node_advection::Union{String, Symbol, Nothing}=nothing,
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    InitializationTools.sync_option_id_with_stype(
        get_options(), model, get_stype_from_model, update_option_id)
    _ = InitializationTools.update_option_id_using_input_option_name(
        get_options(), node_advection, model, get_option_id_from_model, update_option_id)
    dx_topo_api_input = get(kwargs, :dx_topo, nothing)
    dx_topo_current_model = model.topography.parameters.topo_grid.dx_topo.value
    # if dx_topo was loaded from input parameters into the model, update the number of topo nodes
    # Note that loading was either done automatically prior to initialization or
    # via API inputs and the call to load_parameters! above
    update_toponum!(model, dx_topo_current_model)
    # if dx_topo was provided via the API and loaded into the model via load_parameters!, 
    # update the number of topo nodes
    update_toponum!(model, dx_topo_api_input)
    initialize_topography_with_current_parameters!(model)
    update_diff_coeff!(model)
    return nothing
end

function update_toponum!(
    model::ModelData,
    dx_topo::Union{Float64, Nothing}
)::Nothing
    _dx_topo_is_valid = dx_topo_is_valid(dx_topo)
    if !_dx_topo_is_valid
        return nothing
    end
    topo_xsize = model.topography.parameters.topo_grid.topo_xsize.value
    @assert topo_xsize > 0.0 "topo_xsize must be greater than 0.0"
    toponum = calculate_number_of_topo_nodes(topo_xsize, dx_topo)
    @assert toponum > 3 "toponum must be greater than 3"
    print_info("Number of topo nodes: $(toponum) with spacing $(dx_topo) m", level=1)
    model.topography.parameters.topo_grid.toponum.value = toponum
    return nothing
end

function dx_topo_is_valid(dx_topo::Union{Float64, Nothing})::Bool
    return dx_topo !== nothing && dx_topo > 0.0 && !isnan(dx_topo)
end

function initialize_topography_with_current_parameters!(
    model::ModelData
)::Nothing
    update_dx_topo!(model)
    initialize_topo_arrays!(model)
    initialize_topo_and_carb_grids!(model)
    return nothing
end

function get_iuse_topo(model::ModelData)::Int64
    return model.topography.parameters.model_options.iuse_topo.value
end

function get_stype_from_model(model::ModelData)::String
    return model.topography.parameters.model_options.stype_topo.value
end

function get_option_id_from_model(model::ModelData)::Int64
    return model.topography.parameters.model_options.itype_topo.value
end

function update_option_id(model::ModelData, option_id::Int64)::Nothing
    model.topography.parameters.model_options.itype_topo.value = option_id
    return nothing
end

function initialize_topo_and_carb_grids!(
    model::ModelData
)::Nothing
    TopographyFuncs.topography_and_carb_grid_initialize!(model)
    # Deviation from initially flat topography
    TopographyFuncs.topography_and_carb_grid_initialize_triangular_hole!(
        model, check_crustal_hole(model))
    return nothing
end

function check_crustal_hole(model::ModelData)::Bool
    material_description = model.materials.parameters.material_description
    itype_mat = material_description.itype_mat.value
    mat_option_name = OptionTools.get_option_symbol_from_id(mat_option_ids, itype_mat)
    use_crustal_hole = false
    if mat_option_name == mat_option_names.FlexureTriangularHole
        use_crustal_hole = true
    end
    return use_crustal_hole
end

function update_diff_coeff!(
    model::ModelData
)::Nothing
    topo_params = model.topography.parameters
    erosion_velocity = topo_params.depo_and_erosion_rates.erosion_rate.value
    transport_length = topo_params.downhill_diffusion.transport_length.value
    downhill_diff_elev_max = topo_params.downhill_diffusion.downhill_diff_elev_max.value
    topo_diff_coef = erosion_velocity * transport_length^2 / downhill_diff_elev_max
    topo_params.downhill_diffusion.topo_diff_coef.value = topo_diff_coef
    return nothing
end

function update_dx_topo!(
    model::ModelData
)::Nothing
    topo_grid = model.topography.parameters.topo_grid
    toponum = topo_grid.toponum.value
    topo_xsize = topo_grid.topo_xsize.value
    dx_topo = calc_dx_topo(topo_xsize, toponum)
    topo_grid.dx_topo.value = dx_topo
    return nothing
end

""" Check sticky material types if topography is used.

StickyAir and StickyWater material types (mat_type) must be defined if
topography is used. These types are commonly defined in a materials.yml
file but also can be defined interactively.

Note that this method should be called after both topography and marker
materials are initialized.
"""
function check_topography_model(model::ModelData)::Nothing
    if get_iuse_topo(model) == 1
        sticky_air_ids = model.materials.dicts.matid_types["StickyAir"]
        if isempty(sticky_air_ids)
            raise_sticky_error("StickyAir")
        end
        sticky_water_ids = model.materials.dicts.matid_types["StickyWater"]
        if isempty(sticky_water_ids)
            raise_sticky_error("StickyWater")
        end
    end
    return nothing
end

function raise_sticky_error(mat_type::String)::Nothing
    if mat_type ∉ ["StickyWater", "StickyAir"]
        error(
            "Material type $mat_type is not StickyWater or StickyAir in " *
            "sticky error check function."
        )
    end
    error(
        "A material type $mat_type was not found in the materials input " *
        "while using the topography model. Both StickyAir and StickyWater " *
        "material types (mat_type) must be defined when using the topography " *
        "model. Please adjust inputs typically but not always defined in " *
        "materials.yml. Input may also be defined interactively or with a " *
        "materials file with a user defined name."
    )
end

function calc_dx_topo(topo_xsize::Float64, toponum::Int)::Float64
    return topo_xsize / (toponum - 1)
end

function initialize_topo_arrays!(model::ModelData)::Nothing
    toponum = model.topography.parameters.topo_grid.toponum.value
    model.topography.arrays.gridt.array = initialize_topo_array(toponum)
    model.carbonate.arrays.grid_carb.array = initialize_topo_array(toponum)
    return nothing
end

function update_topography_using_velocity_field!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    if get_iuse_topo(model) == 1
        option_id = get_option_id_from_model(model)
        option_name = OptionTools.get_option_symbol_from_id(option_ids, option_id)
        TopoTypeManager.topography_advect!(model, inside_flags, Val(option_name))
    end
    return nothing
end

function calculate_number_of_topo_nodes(xsize::Float64, dx_topo::Float64)::Int64
    result = Int64(floor(xsize/dx_topo)) + 1
    if result < 5
        throw(ArgumentError("The number of topo nodes is less than 5. Decrease dx_topo."))
    end
    return result
end

end # module 