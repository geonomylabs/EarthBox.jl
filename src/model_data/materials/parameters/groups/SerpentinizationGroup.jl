module SerpentinizationGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "serpentinization"

const PDATA = get_eb_parameters()

"""
    Serpentinization <: AbstractParameterGroup

Parameter group for serpentinization model properties.

# Fields
- `iuse_serpentinization::`[`ParameterInt`](@ref): $(PDATA.iuse_serpentinization.description)
- `serpentinization_temperature::`[`ParameterFloat`](@ref): $(PDATA.serpentinization_temperature.description)
- `maximum_serpentinization_depth::`[`ParameterFloat`](@ref): $(PDATA.maximum_serpentinization_depth.description)
- `maximum_serpentinization_rate::`[`ParameterFloat`](@ref): $(PDATA.maximum_serpentinization_rate.description)
- `nominal_strain_rate_serpentinization::`[`ParameterFloat`](@ref): $(PDATA.nominal_strain_rate_serpentinization.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `iuse_serpentinization = $(ROOT_NAME).$(GRP_NAME).iuse_serpentinization.value`
- `serpentinization_temperature = $(ROOT_NAME).$(GRP_NAME).serpentinization_temperature.value`
- `maximum_serpentinization_depth = $(ROOT_NAME).$(GRP_NAME).maximum_serpentinization_depth.value`
- `maximum_serpentinization_rate = $(ROOT_NAME).$(GRP_NAME).maximum_serpentinization_rate.value`
- `nominal_strain_rate_serpentinization = $(ROOT_NAME).$(GRP_NAME).nominal_strain_rate_serpentinization.value`

# Constructor
    Serpentinization()

Create a new Serpentinization parameter group with default values.

# Returns
- `Serpentinization`: New Serpentinization parameter group with initialized values

"""
mutable struct Serpentinization <: AbstractParameterGroup
    iuse_serpentinization::ParameterInt
    serpentinization_temperature::ParameterFloat
    maximum_serpentinization_depth::ParameterFloat
    maximum_serpentinization_rate::ParameterFloat
    nominal_strain_rate_serpentinization::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Serpentinization()::Serpentinization
    pdata = get_eb_parameters()
    data = Serpentinization(
        pdata.iuse_serpentinization,
        pdata.serpentinization_temperature,
        pdata.maximum_serpentinization_depth,
        pdata.maximum_serpentinization_rate,
        pdata.nominal_strain_rate_serpentinization,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
