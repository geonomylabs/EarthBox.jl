module PlotTimeSteppingManager

import ..PlotDtypes: PlotDictType, PlotParametersType

Base.@kwdef mutable struct PlotTimeStepping
    istart::Int = 1
    nsteps::Int = 1
    nskip::Int = 1
    steps::Vector{Int} = Int[]
end

function PlotTimeStepping(plot_dict::PlotDictType)
    istart = plot_dict["general_parameters"]["istart"]
    nsteps = plot_dict["general_parameters"]["nsteps"]
    nskip = plot_dict["general_parameters"]["nskip"]
    steps = get_integer_ids_for_time_steps(istart, nsteps, nskip)
    return PlotTimeStepping(istart=istart, nsteps=nsteps, nskip=nskip, 
                           steps=steps)
end

function set!(
    plot_stepping::PlotTimeStepping, 
    parameters::PlotParametersType
)::Nothing
    if !isnothing(parameters["istart"])
        plot_stepping.istart = parameters["istart"]
    end
    if !isnothing(parameters["nsteps"])
        plot_stepping.nsteps = parameters["nsteps"]
    end
    if !isnothing(parameters["nskip"])
        plot_stepping.nskip = parameters["nskip"]
    end
    plot_stepping.steps = get_integer_ids_for_time_steps(
        plot_stepping.istart, plot_stepping.nsteps, plot_stepping.nskip)
    return nothing
end

function get_integer_ids_for_time_steps(
    istart::Int64, 
    nsteps::Int64, 
    nskip::Int64
)::Vector{Int64}
    steps = Int64[]
    ival = max(istart, 1)
    push!(steps, ival) 
    for i in 1:nsteps
        if i > 1
            ival = ival + nskip + 1
            push!(steps, ival)
        end
    end
    return steps
end

end # module 