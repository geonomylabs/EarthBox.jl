module Parameters

import EarthBox.ModelDataContainer: ModelData
import EarthBox.DictUtils: check_dictionary_names
import EarthBox.ParameterGroupTools: set_group_parameters!
import EarthBox.EarthBoxDtypes: InputDictType
import ..Options: option_names

const key_names = (
    xo_highres = "xo_highres",
    xf_highres = "xf_highres",
    yf_highres = "yf_highres",
    dx_highres = "dx_highres",
    dx_lowres = "dx_lowres",
    dy_highres = "dy_highres",
    dy_lowres = "dy_lowres",
    iuse_trench = "iuse_trench",
    iuse_refinement_delay = "iuse_refinement_delay",
    refinement_time = "refinement_time"
)

const valid_keys = collect(values(key_names))

function set_parameters!(
    model::ModelData, 
    grid_type::Symbol,
    parameters::Union{Dict{String, <:Union{Float64, Int64, String}}, Nothing}
)::Nothing
    if !isnothing(parameters)
        check_dictionary_names(valid_keys, parameters)
        if grid_type == option_names.TtypeRefinedGrid
            set_required_ttype_parameters!(model, parameters)
            set_optional_ttype_parameters!(model, parameters)
        end
    end
    return nothing
end

function set_required_ttype_parameters!(
    model::ModelData,
    parameters::Dict{String, <:Union{Float64, Int64, String}}
)::Nothing
    refinement_param_group = model.grids.parameters.refinement
    # Note that dx_lowres and dy_lowres are not set here because they are not required parameters.
    # These parameters are set in the EarthBoxState constructor and are used to calculate the number 
    # of basic nodes in x and y directions during model initialization.
    parameters_dict = Dict(
        key_names.xo_highres => parameters[key_names.xo_highres],
        key_names.xf_highres => parameters[key_names.xf_highres],
        key_names.dx_highres => parameters[key_names.dx_highres],
        key_names.yf_highres => parameters[key_names.yf_highres],
        key_names.dy_highres => parameters[key_names.dy_highres]
    )
    set_group_parameters!(refinement_param_group, parameters_dict)
    return nothing
end

function set_optional_ttype_parameters!(
    model::ModelData,
    parameters::Dict{String, <:Union{Float64, Int64, String}}
)::Nothing
    refinement_param_group = model.grids.parameters.refinement
    names = [key_names.iuse_trench, key_names.iuse_refinement_delay, key_names.refinement_time]
    for name in names
        if haskey(parameters, name)
            set_group_parameters!(refinement_param_group, Dict(name => parameters[name]))
        end
    end
    return nothing
end

end # module 