module MarkerBoundaryFriction

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    boundary_friction_width::Symbol
    boundary_friction_angle::Symbol
    boundary_cohesion::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize marker boundary friction parameters. 

# Arguments
- `model::ModelData`: The model data container containing the model parameters and arrays.

# Keyword Arguments
- `boundary_friction_width::Union{Float64, Nothing}=nothing`:
    - $(PDATA.boundary_friction_width.description)
- `boundary_friction_angle::Union{Float64, Nothing}=nothing`:
    - $(PDATA.boundary_friction_angle.description)
- `boundary_cohesion::Union{Float64, Nothing}=nothing`:
    - $(PDATA.boundary_cohesion.description)

# Returns
- `Nothing`
"""
function initialize!(
    model::ModelData;
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

end # module