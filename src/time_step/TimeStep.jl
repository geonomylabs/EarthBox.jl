module TimeStep

include("time_step_calculators/TimeStepCalculator.jl")
include("time_step_increase/TimeStepIncrease.jl")
include("update/Update.jl")
include("utils/PrintTimeStep.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit, print_info
import .Update
import .TimeStepIncrease: SingleTimeStepIncrease
import .TimeStepIncrease: MultipleTimeStepIncrease
import .PrintTimeStep

function initialize!(model::ModelData, ntimestep::Int)
    reset_model_time_step_to_viscoelastic_time_step!(model)
    update_time_loop_counter!(model, ntimestep)
    check_and_reset_output_counters!(model)
    print_info(model)
end

function reset_model_time_step_to_viscoelastic_time_step!(model::ModelData)
    main_time_loop = model.timestep.parameters.main_time_loop
    print_info("Resetting main time step to viscoelastic time step", level=2)
    main_time_loop.timestep.value = main_time_loop.timestep_viscoelastic.value
end

function adjust_time_step_using_displacement_limit!(model::ModelData)
    @timeit_memit "Finished updating model time step using displacement limit" begin
        if Update.update_model_time_step_using_displacement_limit(model)
            print_time_step(model, "Updated")
        else
            print_time_step(model, "")
        end
    end
end

function manage_time_step_increase_for_variable_step_models!(model::ModelData)
    ts_params = model.timestep.parameters
    iuse_single = ts_params.single_increase.iuse_single_timestep_increase.value
    iuse_multiple = ts_params.multiple_increase.iuse_multiple_timestep_increase.value
    if iuse_single == 1 || iuse_multiple == 1
        @timeit_memit "Finished managing time step increase for variable step models" begin
            TimeStepIncrease.increase_timestep_and_max_marker_displacement(model)
        end
    end
end

function update_time_loop_counter!(model::ModelData, ntimestep::Int)
    model.timestep.parameters.main_time_loop.ntimestep.value = ntimestep
end

function check_and_reset_output_counters!(model::ModelData)
    icount_output = model.timestep.parameters.output_steps.icount_output
    nskip = model.timestep.parameters.output_steps.nskip
    noutput = model.timestep.parameters.output_steps.noutput

    if icount_output.value == nskip.value
        icount_output.value = 0
    end
    if icount_output.value == 0
        noutput.value += 1
    end
end

function update_output_counter!(model::ModelData)
    model.timestep.parameters.output_steps.icount_output.value += 1
end

function get_maximum_number_of_time_steps(model::ModelData)::Int
    return model.timestep.parameters.main_time_loop.ntimestep_max.value
end

function get_current_number_of_time_steps(model::ModelData)::Int
    return model.timestep.parameters.main_time_loop.ntimestep.value
end

function print_time_step(model::ModelData, msg::String)
    PrintTimeStep.print_timestep(model, msg)
end

function print_info(model::ModelData)
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    sec_per_myr = model.conversion.parameters.sec_per_Myr.value
    timesum_myr = timesum/sec_per_myr
    timestep_myr = timestep/sec_per_myr

    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    icount_output = model.timestep.parameters.output_steps.icount_output.value

    print_info("", level=1)
    print_info(
        "Working on time step: $ntimestep at time (Myr) $(round(timesum_myr; digits=5)) : icount_output: $icount_output : Viscoelastic timestep (Myr): $(round(timestep_myr; digits=5))",
        level=1
    )
    print_info("", level=1)
end

""" Update total model time in seconds (timesum).

Total model time referred to as timesum is advanced using time step updated
based on maximum marker displacement.

Updated Model Parameter
=======================
model.timestep.parameters.main_time_loop
----------------------------------------
timesum.value: float
    - Total model time in seconds.
"""
function advance_model_time!(model::ModelData)::Nothing
    timestep = model.timestep.parameters.main_time_loop.timestep
    timesum = model.timestep.parameters.main_time_loop.timesum
    timesum.value = timesum.value + timestep.value
    return nothing
end

""" Print initial and current sea level coordinates.

Parameters
----------
model::ModelData
    - The model data containing sea level information.
"""
function print_initial_sealevel(model::ModelData)::Nothing
    sealevel = model.topography.parameters.sealevel
    y_water_ini = sealevel.y_water_ini.value
    y_sealevel = sealevel.y_sealevel.value
    print_info("Initial water level y-coordinate (m) : $y_water_ini", level=1)
    print_info("Current sealevel y-coordinate (m) : $y_sealevel", level=1)
    return nothing
end

""" Print loop information including completion time for each time step.

Parameters
----------
ntimestep::Int
    - The current time step number.
loop_time_start::Float64
    - The start time of the loop.
"""
function print_loop_info(ntimestep::Int, loop_time_start::Float64)::Nothing
    time2 = time()
    print_info("Completed time step $ntimestep: cpu(s) $(round(time2 - loop_time_start, digits=4))", level=1)
    print_info("", level=1)
    return nothing
end

""" Set time loop parameters.

Parameters
----------
model::ModelData
    - The model data.
number_of_time_steps::Union{Int, Nothing}
    - The number of time steps.
viscoelastic_time_step_seconds::Union{Float64, Nothing}
    - The viscoelastic time step in seconds.
output_time_step_seconds::Union{Float64, Nothing}
    - The output time step in seconds.
"""
function set_parameters!(
    model::ModelData;
    number_of_time_steps::Union{Int, Nothing}=nothing,
    viscoelastic_time_step_seconds::Union{Float64, Nothing}=nothing,
    output_time_step_seconds::Union{Float64, Nothing}=nothing
)::Nothing
    main_time_loop = model.timestep.parameters.main_time_loop
    if !isnothing(number_of_time_steps)
        main_time_loop.ntimestep_max.value = number_of_time_steps
    end
    if !isnothing(viscoelastic_time_step_seconds)
        main_time_loop.timestep_viscoelastic.value = viscoelastic_time_step_seconds
    end
    output_steps = model.timestep.parameters.output_steps
    if !isnothing(output_time_step_seconds)
        output_steps.timestep_out.value = output_time_step_seconds
    end
    return nothing
end

""" Update nskip parameter in output_steps.

Parameters
----------
model::ModelData
    - The model data.
"""
function update_nskip!(model::ModelData)::Nothing
    timestep_params = model.timestep.parameters
    timestep_out = timestep_params.output_steps.timestep_out.value
    timestep_viscoelastic = timestep_params.main_time_loop.timestep_viscoelastic.value
    nskip = floor(Int, timestep_out/timestep_viscoelastic)
    timestep_params.output_steps.nskip.value = nskip
    return nothing
end

""" Update model time step.

Parameters
----------
model::ModelData
    - The model data.
"""
function update_model_time_step!(model::ModelData)::Nothing
    main_time_loop = model.timestep.parameters.main_time_loop
    timestep_viscoelastic = main_time_loop.timestep_viscoelastic.value
    main_time_loop.timestep.value = timestep_viscoelastic
    return nothing
end

end # module 