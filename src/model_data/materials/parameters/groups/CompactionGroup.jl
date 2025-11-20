module CompactionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.materials.parameters"
const GRP_NAME = "compaction"

const PDATA = get_eb_parameters()

"""
    Compaction <: AbstractParameterGroup

Parameter group for sediment compaction and water properties.

# Fields
- `iuse_sed_porosity::`[`ParameterInt`](@ref): $(PDATA.iuse_sed_porosity.description)
- `conductivity_water::`[`ParameterFloat`](@ref): $(PDATA.conductivity_water.description)
- `density_water::`[`ParameterFloat`](@ref): $(PDATA.density_water.description)
- `heat_capacity_water::`[`ParameterFloat`](@ref): $(PDATA.heat_capacity_water.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of 
parameter objects

# Nested Dot Access
- `iuse_sed_porosity = $(ROOT_NAME).$(GRP_NAME).iuse_sed_porosity.value`
- `conductivity_water = $(ROOT_NAME).$(GRP_NAME).conductivity_water.value`
- `density_water = $(ROOT_NAME).$(GRP_NAME).density_water.value`
- `heat_capacity_water = $(ROOT_NAME).$(GRP_NAME).heat_capacity_water.value`

# Constructor
    Compaction()

Create a new Compaction parameter group with default values.

# Returns
- `Compaction`: New Compaction parameter group with initialized values

"""
mutable struct Compaction <: AbstractParameterGroup
    iuse_sed_porosity::ParameterInt
    conductivity_water::ParameterFloat
    density_water::ParameterFloat
    heat_capacity_water::ParameterFloat
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Compaction()::Compaction
    pdata = get_eb_parameters()
    data = Compaction(
        pdata.iuse_sed_porosity,
        pdata.conductivity_water,
        pdata.density_water,
        pdata.heat_capacity_water,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
