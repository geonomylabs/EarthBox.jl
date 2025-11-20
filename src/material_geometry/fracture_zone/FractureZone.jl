module FractureZone

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!

const PDATA = get_eb_parameters()

struct ValidInputNames
    sediment_thickness::Symbol
    basaltic_oceanic_crust_thickness::Symbol
    gabbroic_oceanic_crust_thickness::Symbol
    thickness_of_younger_lithosphere::Symbol
    thickness_of_older_lithosphere::Symbol
    thickness_of_weak_lithosphere::Symbol
    x_fracture_zone_start::Symbol
    x_fracture_zone_end::Symbol
end

""" Initialize fracture zone geometry.

# Keyword arguments:
- `sediment_thickness::Float64`:
    - $(PDATA.sediment_thickness.description)
- `basaltic_oceanic_crust_thickness::Float64`:
    - $(PDATA.basaltic_oceanic_crust_thickness.description)
- `gabbroic_oceanic_crust_thickness::Float64`:
    - $(PDATA.gabbroic_oceanic_crust_thickness.description)
- `thickness_of_younger_lithosphere::Float64`:
    - $(PDATA.thickness_of_younger_lithosphere.description)
- `thickness_of_older_lithosphere::Float64`:
    - $(PDATA.thickness_of_older_lithosphere.description)
- `thickness_of_weak_lithosphere::Float64`:
    - $(PDATA.thickness_of_weak_lithosphere.description)
- `x_fracture_zone_start::Float64`:
    - $(PDATA.x_fracture_zone_start.description)
- `x_fracture_zone_end::Float64`:
    - $(PDATA.x_fracture_zone_end.description)
"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

end 