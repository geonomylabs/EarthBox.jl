module Pressure

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ModelDataContainer: set_parameter!

const PDATA = get_eb_parameters()

struct ValidInputNames
    pressure_bc::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize pressure boundary conditions. The pressure boundary condition is either applied in the upper 
left corner of the model or along the top or bottom boundaries depending on the boundary
condition type selected by the user. See [`BoundaryConditions.initialize!`](@ref EarthBox.BoundaryConditions.initialize!) 
for more information.

# Arguments
- `model::`[`ModelData`](@ref ModelData): Model data object containing model parameters and arrays.

# Optional Keyword Arguments
- `pressure_bc::Union{Float64, Nothing}=nothing`: 
    - $(PDATA.pressure_bc.description)

# Returns
- `Nothing`

"""
function initialize!(model::ModelData; kwargs...)
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

end # module