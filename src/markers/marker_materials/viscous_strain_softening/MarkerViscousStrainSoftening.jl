module MarkerViscousStrainSoftening

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.ModelDataContainer: set_parameter!
import EarthBox.UseOptionTools: set_use_option!

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_viscous_strain_soft::Symbol
    vsoftfac::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize marker viscous strain softening for plasticity models. This initializes 
the `iuse_viscous_strain_soft` and `vsoftfac` parameters in the model data container 
and is required for all models that use viscous strain softening whereby the pre-exponential
term for dislocation creep is modified by a function of plastic strain and the strain softening 
factor `vsoftfac`.

# Arguments
- `model::ModelData`: The model data container containing the model parameters and arrays.

# Keyword Arguments
- `iuse_viscous_strain_soft::Union{Int, Nothing}=nothing`:
    - $(PDATA.iuse_viscous_strain_soft.description)
- `vsoftfac::Union{Float64, Nothing}=nothing`:
    - $(PDATA.vsoftfac.description)

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