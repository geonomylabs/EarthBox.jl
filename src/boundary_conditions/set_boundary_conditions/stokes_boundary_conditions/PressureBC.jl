"""
    Function used to set pressure bc mode.
"""
module PressureBC

import EarthBox.ModelDataContainer: ModelData

"""
    set_pressure_bc_mode!(model::ModelData, mode::Int)::Nothing

Set the pressure boundary condition mode.

# Arguments
- `model::ModelData`: Model data.
- `mode::Int`: Pressure boundary condition mode. Options are:
    - 0: pressure prescribed at upper left node.
    - 1: pressure prescribed along top and bottom boundaries with top pressure set
         to the prescribed pressure value and the bottom pressure set to zero.
    - 2: pressure prescribed along left and right boundaries with the left pressure
         set to the prescribed pressure value and the right pressure set to zero.

The application of the pressure boundary condition is done during the assembly of
the linear system of equations.
"""
function set_pressure_bc_mode!(model::ModelData, mode::Int)::Nothing
    model.bcs.parameters.bc_options.pressure_bc_mode.value = mode
    return nothing
end

end # module
