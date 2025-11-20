module MarkerAndGridInterpOptionsGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.interpolation.parameters"
const GRP_NAME = "interp_options"

const PDATA = get_eb_parameters()

"""
    MarkerAndGridInterpOptions <: AbstractParameterGroup

Parameter group for marker-grid interpolation options.

# Fields
- `iuse_initial_temp_interp::`[`ParameterInt`](@ref): $(PDATA.iuse_initial_temp_interp.description)
- `iuse_harmonic_avg_normal_viscosity::`[`ParameterInt`](@ref): $(PDATA.iuse_harmonic_avg_normal_viscosity.description)
- `obj_list::Vector{ParameterInt}`: List of parameter objects

# Nested Dot Access
- `iuse_initial_temp_interp = $(ROOT_NAME).$(GRP_NAME).iuse_initial_temp_interp.value`
- `iuse_harmonic_avg_normal_viscosity = $(ROOT_NAME).$(GRP_NAME).iuse_harmonic_avg_normal_viscosity.value`

# Constructor
    MarkerAndGridInterpOptions()

Create a new MarkerAndGridInterpOptions parameter group with default values.

# Returns
- `MarkerAndGridInterpOptions`: New parameter group with initialized values

"""
mutable struct MarkerAndGridInterpOptions <: AbstractParameterGroup
    iuse_initial_temp_interp::ParameterInt
    iuse_harmonic_avg_normal_viscosity::ParameterInt
    obj_list::Vector{ParameterInt}
end

function MarkerAndGridInterpOptions()::MarkerAndGridInterpOptions
    pdata = get_eb_parameters()
    data = MarkerAndGridInterpOptions(
        pdata.iuse_initial_temp_interp,
        pdata.iuse_harmonic_avg_normal_viscosity,
        ParameterInt[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
