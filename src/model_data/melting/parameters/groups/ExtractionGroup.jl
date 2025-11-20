module ExtractionGroup

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.Parameters: ParameterFloat, ParameterInt
import EarthBox.Parameters: print_parameter
import EarthBox.ParameterGroupTools: get_numerical_parameter_object_list
import EarthBox.EarthBoxDtypes: AbstractParameterGroup

const ROOT_NAME = "model.melting.parameters"
const GRP_NAME = "extraction"

const PDATA = get_eb_parameters()

"""
    Extraction <: AbstractParameterGroup

Parameter group for melt extraction model parameters.

# Fields
- `melt_residual::`[`ParameterFloat`](@ref): $(PDATA.melt_residual.description)
- `ext_vol::`[`ParameterFloat`](@ref): $(PDATA.ext_vol.description)
- `xmid_mol::`[`ParameterFloat`](@ref): $(PDATA.xmid_mol.description)
- `ytop_mol::`[`ParameterFloat`](@ref): $(PDATA.ytop_mol.description)
- `width_mol::`[`ParameterFloat`](@ref): $(PDATA.width_mol.description)
- `ndrainage_basin::`[`ParameterInt`](@ref): $(PDATA.ndrainage_basin.description)
- `smoothing_radius_drainage::`[`ParameterFloat`](@ref): $(PDATA.smoothing_radius_drainage.description)
- `characteristic_injection_width::`[`ParameterFloat`](@ref): $(PDATA.characteristic_injection_width.description)
- `magma_height_limit::`[`ParameterFloat`](@ref): $(PDATA.magma_height_limit.description)
- `fractionation_threshold_limit::`[`ParameterFloat`](@ref): $(PDATA.fractionation_threshold_limit.description)
- `emplacement_temperature::`[`ParameterFloat`](@ref): $(PDATA.emplacement_temperature.description)
- `number_of_injection_subdomains::`[`ParameterInt`](@ref): $(PDATA.number_of_injection_subdomains.description)
- `maximum_shallow_injection_depth::`[`ParameterFloat`](@ref): $(PDATA.maximum_shallow_injection_depth.description)
- `extraction_fraction::`[`ParameterFloat`](@ref): $(PDATA.extraction_fraction.description)
- `smoothing_radius_fractionation::`[`ParameterFloat`](@ref): $(PDATA.smoothing_radius_fractionation.description)
- `mantle_search_width::`[`ParameterFloat`](@ref): $(PDATA.mantle_search_width.description)
- `ndrainage_basin_old::`[`ParameterInt`](@ref): $(PDATA.ndrainage_basin_old.description)
- `obj_list::Vector{Union{ParameterFloat, ParameterInt}}`: List of parameter objects

# Nested Dot Access
- `melt_residual = $(ROOT_NAME).$(GRP_NAME).melt_residual.value`
- `ext_vol = $(ROOT_NAME).$(GRP_NAME).ext_vol.value`
- `xmid_mol = $(ROOT_NAME).$(GRP_NAME).xmid_mol.value`
- `ytop_mol = $(ROOT_NAME).$(GRP_NAME).ytop_mol.value`
- `width_mol = $(ROOT_NAME).$(GRP_NAME).width_mol.value`
- `ndrainage_basin = $(ROOT_NAME).$(GRP_NAME).ndrainage_basin.value`
- `smoothing_radius_drainage = $(ROOT_NAME).$(GRP_NAME).smoothing_radius_drainage.value`
- `characteristic_injection_width = $(ROOT_NAME).$(GRP_NAME).characteristic_injection_width.value`
- `magma_height_limit = $(ROOT_NAME).$(GRP_NAME).magma_height_limit.value`
- `fractionation_threshold_limit = $(ROOT_NAME).$(GRP_NAME).fractionation_threshold_limit.value`
- `emplacement_temperature = $(ROOT_NAME).$(GRP_NAME).emplacement_temperature.value`
- `number_of_injection_subdomains = $(ROOT_NAME).$(GRP_NAME).number_of_injection_subdomains.value`
- `maximum_shallow_injection_depth = $(ROOT_NAME).$(GRP_NAME).maximum_shallow_injection_depth.value`
- `extraction_fraction = $(ROOT_NAME).$(GRP_NAME).extraction_fraction.value`
- `smoothing_radius_fractionation = $(ROOT_NAME).$(GRP_NAME).smoothing_radius_fractionation.value`
- `mantle_search_width = $(ROOT_NAME).$(GRP_NAME).mantle_search_width.value`
- `ndrainage_basin_old = $(ROOT_NAME).$(GRP_NAME).ndrainage_basin_old.value`

# Constructor
    Extraction()

# Returns
- `Extraction`: New Extraction parameter group with initialized values

"""
mutable struct Extraction <: AbstractParameterGroup
    melt_residual::ParameterFloat
    ext_vol::ParameterFloat
    xmid_mol::ParameterFloat
    ytop_mol::ParameterFloat
    width_mol::ParameterFloat
    ndrainage_basin::ParameterInt
    smoothing_radius_drainage::ParameterFloat
    characteristic_injection_width::ParameterFloat
    magma_height_limit::ParameterFloat
    fractionation_threshold_limit::ParameterFloat
    emplacement_temperature::ParameterFloat
    number_of_injection_subdomains::ParameterInt
    maximum_shallow_injection_depth::ParameterFloat
    extraction_fraction::ParameterFloat
    smoothing_radius_fractionation::ParameterFloat
    mantle_search_width::ParameterFloat
    ndrainage_basin_old::ParameterInt
    obj_list::Vector{Union{ParameterFloat, ParameterInt}}
end

function Extraction()::Extraction
    pdata = get_eb_parameters()
    data = Extraction(
        pdata.melt_residual,
        pdata.ext_vol,
        pdata.xmid_mol,
        pdata.ytop_mol,
        pdata.width_mol,
        pdata.ndrainage_basin,
        pdata.smoothing_radius_drainage,
        pdata.characteristic_injection_width,
        pdata.magma_height_limit,
        pdata.fractionation_threshold_limit,
        pdata.emplacement_temperature,
        pdata.number_of_injection_subdomains,
        pdata.maximum_shallow_injection_depth,
        pdata.extraction_fraction,
        pdata.smoothing_radius_fractionation,
        pdata.mantle_search_width,
        pdata.ndrainage_basin_old,
        Union{ParameterFloat, ParameterInt}[] # obj_list
    )
    data.obj_list = get_numerical_parameter_object_list(data)
    return data
end

end # module
