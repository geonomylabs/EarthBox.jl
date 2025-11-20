module TimeLoopUtils

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData

function update_main_time_counters!(model::ModelData, _ntimestep::Int)::Nothing
    model.timestep.parameters.main_time_loop.ntimestep.value = _ntimestep
    return nothing
end

function print_timestep_info(model::ModelData)::Nothing
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    sec_per_myr = model.conversion.parameters.sec_per_Myr.value
    timesum_myr = timesum/sec_per_myr
    timestep_myr = timestep/sec_per_myr

    _ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    icount_output = model.timestep.parameters.output_steps.icount_output.value

    print_info("", level=1)
    print_info(
        "Working on time step: $_ntimestep at time (Myr) $(round(timesum_myr; digits=5)) : "
        *"icount_output: $icount_output : Viscoelastic timestep (Myr): $(round(timestep_myr; digits=5))",
        level=1
    )
    return nothing
end

function check_and_reset_output_counters!(model::ModelData)::Nothing
    icount_output = model.timestep.parameters.output_steps.icount_output
    nskip = model.timestep.parameters.output_steps.nskip
    noutput = model.timestep.parameters.output_steps.noutput

    if icount_output.value == nskip.value
        icount_output.value = 0
    end
    if icount_output.value == 0
        noutput.value += 1
    end
    return nothing
end

function update_output_counter!(model::ModelData)::Nothing
    model.timestep.parameters.output_steps.icount_output.value += 1
    return nothing
end

end # module