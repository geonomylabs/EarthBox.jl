module TimeStepIncrease

include("types/MultipleTimeStepIncrease.jl")
include("types/SingleTimeStepIncrease.jl")

import EarthBox.ModelDataContainer: ModelData

function increase_timestep_and_max_marker_displacement(model::ModelData)
    ntime_increase = model.timestep.parameters.single_increase.ntime_increase.value

    ts_params = model.timestep.parameters
    ntimestep = ts_params.main_time_loop.ntimestep.value

    iuse_single_timestep_increase = 
        ts_params.single_increase.iuse_single_timestep_increase.value
    iuse_multiple_timestep_increase = 
        ts_params.multiple_increase.iuse_multiple_timestep_increase.value

    if iuse_single_timestep_increase == 1
        # Only update timestep for single increase model
        if ntimestep == ntime_increase
            update_timestep(model)
        end
    elseif iuse_multiple_timestep_increase == 1
        # Update timestep and max displacement for multiple increase model
        if ntimestep in get_change_points(model)
            update_timestep(model)
            update_max_displacement(model)
            print_info(model)
        end
    end
end

"""
    get_change_points(model::ModelData)::Vector{Int}

Get change points for multiple increase model.

Get time steps where time step and max marker displacement are increased.
"""
function get_change_points(model::ModelData)::Vector{Int}
    ts_params = model.timestep.parameters
    ntime_increase_1 = ts_params.multiple_increase.ntime_increase_1.value
    ntime_increase_2 = ts_params.multiple_increase.ntime_increase_2.value
    ntime_increase_3 = ts_params.multiple_increase.ntime_increase_3.value
    ntime_increase_4 = ts_params.multiple_increase.ntime_increase_4.value
    return [ntime_increase_1, ntime_increase_2, ntime_increase_3, ntime_increase_4]
end

function print_info(model::ModelData)
    ts_params = model.timestep.parameters
    timestep_viscoelastic = ts_params.main_time_loop.timestep_viscoelastic
    marker_cell_displ_max = 
        model.markers.parameters.advection.marker_cell_displ_max
    println("\n >>>>>> timestep_viscoelastic, marker_cell_displ_max : ",
            timestep_viscoelastic.value, " ", marker_cell_displ_max.value)
end

function update_timestep(model::ModelData)
    factor = model.timestep.parameters.multiple_increase.time_increase_factor.value
    ts_params = model.timestep.parameters
    timestep_viscoelastic = ts_params.main_time_loop.timestep_viscoelastic
    timestep_viscoelastic.value = timestep_viscoelastic.value * factor
end

function update_max_displacement(model::ModelData)
    factor = model.timestep.parameters.multiple_increase.cell_displ_factor.value
    marker_cell_displ_max = model.markers.parameters.advection.marker_cell_displ_max
    marker_cell_displ_max.value = marker_cell_displ_max.value * factor
end

end # module 