module RockProperties

include("adiabatic_term/UpdateAdiabaticTerm.jl")
include("density/DensityModel.jl")
include("rhocp/RhoCpModel.jl")
include("thermal_conductivity/ThermalConductivityModel.jl")
include("porosity/PorosityUpdate.jl")

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.EarthBoxDtypes: InputDictType
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.UseOptionTools: set_use_option!
import EarthBox.ParameterGroupTools: set_group_parameters!
import EarthBox: InitializationTools
import EarthBox: OptionTools
import .DensityModel
import .RhoCpModel
import .ThermalConductivityModel
import .UpdateAdiabaticTerm
import .PorosityUpdate

export initialize!

const PDATA = get_eb_parameters()

struct ValidInputNames
    maximum_heat_capacity::Symbol
    iuse_sed_porosity::Symbol
    iuse_melt_lens::Symbol
end

"""
    initialize!(model::ModelData)::Nothing

Initialize rock properties for all markers.

# Arguments
- `model::ModelData`: Model data container containing the model parameters and arrays.

- `density_model::Union{String, Symbol, Nothing}=nothing`: 
    Density model option name. See the **Density Models** section below for 
    information on available density models.

- `thermal_conductivity_model::Union{String, Symbol, Nothing}=nothing`: 
    Thermal conductivity model option name. See the **Thermal Conductivity Models** 
    section below for information on available thermal conductivity models.

- `rhocp_model::Union{String, Symbol, Nothing}=nothing`: 
    Heat capacity model option name. See the **RhoCp Models** section below for 
    information on available heat capacity models.

# Keyword arguments
- `iuse_sed_porosity::Int64`: 
    - $(PDATA.iuse_sed_porosity.description)

- `iuse_melt_lens::Int64`: 
    - $(PDATA.iuse_melt_lens.description)

- `maximum_heat_capacity::Float64`: 
    - $(PDATA.maximum_heat_capacity.description)
---
# Density Models
---
$(DensityModel.make_density_models_string())
---
# Thermal Conductivity Models
---
$(ThermalConductivityModel.make_thermal_conductivity_models_string())
---
# RhoCp Models
---
$(RhoCpModel.make_rhocp_models_string())

"""
function initialize!(
    model::ModelData;
    density_model::Union{String, Symbol, Nothing}=nothing,
    thermal_conductivity_model::Union{String, Symbol, Nothing}=nothing,
    rhocp_model::Union{String, Symbol, Nothing}=nothing,
    kwargs...
)::Nothing
    DensityModel.initialize!(
        model, density_model=density_model)
    ThermalConductivityModel.initialize!(
        model, thermal_conductivity_model=thermal_conductivity_model)
    RhoCpModel.initialize!(
        model, rhocp_model=rhocp_model)
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

function print_option(model::ModelData)::Nothing
    DensityModel.print_option(model)
    ThermalConductivityModel.print_option(model)
    RhoCpModel.print_option(model)
    return nothing
end

function update_marker_rock_properties!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    DensityModel.update_density!(model, inside_flags)
    RhoCpModel.update_rhocp!(model, inside_flags)
    ThermalConductivityModel.update_conductivity!(model, inside_flags)
    @timeit_memit "Finished updating marker adiabatic term" begin
        UpdateAdiabaticTerm.update_adiabatic_term!(model, inside_flags)
    end
    update_for_porosity!(model, inside_flags)
    return nothing
end

function update_for_porosity!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    iuse_sed_porosity = model.materials.parameters.compaction.iuse_sed_porosity.value
    if iuse_sed_porosity == 1
        @timeit_memit "Finished updating marker porosity" begin
            PorosityUpdate.update_rock_properties_for_porosity!(model, inside_flags)
        end
    end
    return nothing
end

end # module 