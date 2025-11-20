module MeltDamageGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "melt_damage"

const PDATA = get_eb_parameters()

"""
    MeltDamage <: AbstractParameterGroup

Parameter group for melt damage properties.

# Fields
- `iuse_melt_damage::`[`ParameterInt`](@ref): $(PDATA.iuse_melt_damage.description)
- `melt_damage_distance::`[`ParameterFloat`](@ref): $(PDATA.melt_damage_distance.description)
- `melt_damage_factor::`[`ParameterFloat`](@ref): $(PDATA.melt_damage_factor.description)
- `melt_damage_taper_distance::`[`ParameterFloat`](@ref): $(PDATA.melt_damage_taper_distance.description)
- `iuse_probabilistic_melt_damage::`[`ParameterInt`](@ref): $(PDATA.iuse_probabilistic_melt_damage.description)
- `maximum_damage_probability::`[`ParameterFloat`](@ref): $(PDATA.maximum_damage_probability.description)
- `magmatic_crust_height_threshold::`[`ParameterFloat`](@ref): $(PDATA.magmatic_crust_height_threshold.description)
- `magmatic_crust_height_minimum::`[`ParameterFloat`](@ref): $(PDATA.magmatic_crust_height_minimum.description)
- `magmatic_crust_height_maximum::`[`ParameterFloat`](@ref): $(PDATA.magmatic_crust_height_maximum.description)
- `magmatic_crust_height_intermediate::`[`ParameterFloat`](@ref): $(PDATA.magmatic_crust_height_intermediate.description)
- `intermediate_damage_probability::`[`ParameterFloat`](@ref): $(PDATA.intermediate_damage_probability.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of 
parameter objects

# Nested Dot Access
- `iuse_melt_damage = $(ROOT_NAME).$(GRP_NAME).iuse_melt_damage.value`
- `melt_damage_distance = $(ROOT_NAME).$(GRP_NAME).melt_damage_distance.value`
- `melt_damage_factor = $(ROOT_NAME).$(GRP_NAME).melt_damage_factor.value`
- `melt_damage_taper_distance = $(ROOT_NAME).$(GRP_NAME).melt_damage_taper_distance.value`
- `iuse_probabilistic_melt_damage = $(ROOT_NAME).$(GRP_NAME).iuse_probabilistic_melt_damage.value`
- `maximum_damage_probability = $(ROOT_NAME).$(GRP_NAME).maximum_damage_probability.value`
- `magmatic_crust_height_threshold = $(ROOT_NAME).$(GRP_NAME).magmatic_crust_height_threshold.value`
- `magmatic_crust_height_minimum = $(ROOT_NAME).$(GRP_NAME).magmatic_crust_height_minimum.value`
- `magmatic_crust_height_maximum = $(ROOT_NAME).$(GRP_NAME).magmatic_crust_height_maximum.value`
- `magmatic_crust_height_intermediate = $(ROOT_NAME).$(GRP_NAME).magmatic_crust_height_intermediate.value`
- `intermediate_damage_probability = $(ROOT_NAME).$(GRP_NAME).intermediate_damage_probability.value`

# Constructor
    MeltDamage()

# Returns
- `MeltDamage`: New MeltDamage parameter group with initialized values

"""
mutable struct MeltDamage <: AbstractParameterGroup
    iuse_melt_damage::ParameterInt
    melt_damage_distance::ParameterFloat
    melt_damage_factor::ParameterFloat
    melt_damage_taper_distance::ParameterFloat
    iuse_probabilistic_melt_damage::ParameterInt
    maximum_damage_probability::ParameterFloat
    magmatic_crust_height_threshold::ParameterFloat
    magmatic_crust_height_minimum::ParameterFloat
    magmatic_crust_height_maximum::ParameterFloat
    magmatic_crust_height_intermediate::ParameterFloat
    intermediate_damage_probability::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function MeltDamage()::MeltDamage
    pdata = get_eb_parameters()
    data = MeltDamage(
        pdata.iuse_melt_damage,
        pdata.melt_damage_distance,
        pdata.melt_damage_factor,
        pdata.melt_damage_taper_distance,
        pdata.iuse_probabilistic_melt_damage,
        pdata.maximum_damage_probability,
        pdata.magmatic_crust_height_threshold,
        pdata.magmatic_crust_height_minimum,
        pdata.magmatic_crust_height_maximum,
        pdata.magmatic_crust_height_intermediate,
        pdata.intermediate_damage_probability,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
