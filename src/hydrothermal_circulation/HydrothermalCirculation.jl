module HydrothermalCirculation

include("core/UpdateHydrothermalCirculation.jl")

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ConversionFuncs: kelvin_to_celsius

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_hydrothermal::Symbol
    iuse_plastic_strain_rate_for_hydrothermal::Symbol
    iuse_plastic_strain_for_hydrothermal::Symbol
    hydrothermal_smoothing_factor::Symbol
    hydrothermal_nusselt_number::Symbol
    hydrothermal_max_temperature::Symbol
    hydrothermal_max_depth::Symbol
    hydrothermal_decay_length::Symbol
    hydrothermal_buffer_distance::Symbol
    hydrothermal_plastic_strain_rate_reference::Symbol
    hydrothermal_plastic_strain_reference::Symbol
    sediment_thickness_threshold::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize parameters for the hydrothermal circulation model used to approximate
the effects of hydrothermal circulation on thermal structure using an effective
thermal conductivity.

# Keyword Arguments
- `$(PDATA.iuse_hydrothermal.name)::Int64`
    - $(PDATA.iuse_hydrothermal.description)
- `$(PDATA.iuse_plastic_strain_rate_for_hydrothermal.name)::Int64`
    - $(PDATA.iuse_plastic_strain_rate_for_hydrothermal.description)
- `$(PDATA.iuse_plastic_strain_for_hydrothermal.name)::Int64`
    - $(PDATA.iuse_plastic_strain_for_hydrothermal.description)
- `$(PDATA.hydrothermal_smoothing_factor.name)::Float64`
    - $(PDATA.hydrothermal_smoothing_factor.description)
- `$(PDATA.hydrothermal_nusselt_number.name)::Float64`
    - $(PDATA.hydrothermal_nusselt_number.description)
- `$(PDATA.hydrothermal_max_temperature.name)::Float64`
    - $(PDATA.hydrothermal_max_temperature.description)
- `$(PDATA.hydrothermal_max_depth.name)::Float64`
    - $(PDATA.hydrothermal_max_depth.description)
- `$(PDATA.hydrothermal_decay_length.name)::Float64`
    - $(PDATA.hydrothermal_decay_length.description)
- `$(PDATA.hydrothermal_buffer_distance.name)::Float64`
    - $(PDATA.hydrothermal_buffer_distance.description)
- `$(PDATA.hydrothermal_plastic_strain_rate_reference.name)::Float64`
    - $(PDATA.hydrothermal_plastic_strain_rate_reference.description)
- `$(PDATA.hydrothermal_plastic_strain_reference.name)::Float64`
    - $(PDATA.hydrothermal_plastic_strain_reference.description)
- `$(PDATA.sediment_thickness_threshold.name)::Float64`
    - $(PDATA.sediment_thickness_threshold.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

function update_rock_properties_for_hydrothermal!(
    model::ModelData
)::Nothing
    iuse_hydrothermal = model.materials.parameters.hydrothermal.iuse_hydrothermal.value
    if iuse_hydrothermal == 1
        @timeit_memit "Finished updating marker rock properties for hydrothermal circulation" begin
            UpdateHydrothermalCirculation.update_for_hydrothermal_opt!(model)
        end
    end
    return nothing
end


end # module