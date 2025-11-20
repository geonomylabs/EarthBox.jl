module WeakFault

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    fault_dip_degrees::Symbol
    fault_thickness::Symbol
    x_initial_fault::Symbol
    fault_height::Symbol
    iuse_weak_fault::Symbol
end

""" Initialize weak fault geometry.

# Keyword arguments:
- `fault_dip_degrees::Float64`:
    - $(PDATA.fault_dip_degrees.description)
- `fault_thickness::Float64`:
    - $(PDATA.fault_thickness.description)
- `x_initial_fault::Float64`:
    - $(PDATA.x_initial_fault.description)
- `fault_height::Float64`:
    - $(PDATA.fault_height.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    # Turn on weak fault when initialization geometry is used
    model.geometry.parameters.weak_fault.iuse_weak_fault.value = 1
    return nothing
end

end 