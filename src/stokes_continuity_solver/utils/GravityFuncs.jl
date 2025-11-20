module GravityFuncs

import Printf
import EarthBox.ModelDataContainer: ModelData

"""
    turn_off_gravity(model::ModelData)

Set the vertical component of gravity to zero.
"""
function turn_off_gravity!(model::ModelData)
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    grav_params = model.gravity.parameters
    nsteps_turn_off_gravity = grav_params.nsteps_turn_off_gravity.value
    
    gravity_y = grav_params.gravity_y
    if ntimestep == nsteps_turn_off_gravity
        gravity_y.value = 0.0
    end
end

end # module 