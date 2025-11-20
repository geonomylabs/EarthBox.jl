module UseOptionTools

import EarthBox.EarthBoxDtypes: AbstractParameterGroup
import EarthBox.EarthBoxDtypes: InputDictType
import EarthBox.ParameterGroupTools: get_parameter_value
import EarthBox.ParameterGroupTools: set_group_parameters!

function set_use_option!(
    use_option::Union{Bool, Nothing},
    parameter_group::AbstractParameterGroup,
    iuse_name::String
)::Union{Bool, Nothing}
    iuse_init = get_parameter_value(parameter_group, iuse_name)
    if use_option === nothing
        use_option = get_use_bool_from_iuse(iuse_init)
    else
        iuse = get_iuse_from_use_bool(use_option)
        set_group_parameters!(parameter_group, InputDictType(iuse_name => iuse))
    end
    return use_option
end

function get_use_bool_from_iuse(iuse::Int)::Bool
    return iuse == 1
end

function get_iuse_from_use_bool(use_bool::Bool)::Int
    return use_bool ? 1 : 0
end

end # module UseOptionTools 