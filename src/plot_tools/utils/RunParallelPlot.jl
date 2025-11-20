module RunParallelPlot

export run_parallel_marker_plotter

"""
    run_parallel_marker_plotter()::Nothing

Run the marker plotting script `Plot.jl` in parallel for a given range of time steps.

The script `Plot.jl` must be located in the model directory along with the script
that executes this function. The script `Plot.jl` is called via command line as 
follows:

```bash
julia Plot.jl marker_plots case_name=<case_name> istart=<istart> iend=<iend>
```

This assumes that the command line plotter function `run_cl_plotter()` is being
passed the root path to the model output directory when called from the `Plot.jl` 
script.

Usage:

```bash`
julia `RunParallelPlot.jl` <case_name> <total_steps> <num_processors>
```

where:
- <case_name>: The name of the case.
- <total_steps>: The total number of time steps.
- <num_processors>: The number of processors to use.

and `RunParallelPlot.jl` is the script that imports and executes this function
and is located in the model directory along with the `Plot.jl` script.
An example of the `RunParallelPlot.jl` script is as follows:

```julia
module RunParallelPlot

using EarthBox

function main()::Nothing
    run_parallel_marker_plotter()
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    RunParallelPlot.main()
end
```
"""
function run_parallel_marker_plotter()::Nothing
    if length(ARGS) != 3
        error("There were too few arguments. Usage: julia ParallelMarkerPlotter.jl <case_name> <total_steps> <num_processors>")
    end
    case_name = ARGS[1]
    if occursin("case", case_name) == false
        error("case_name must contain the substring \"case\" (e.g., \"case0\", \"case1\").")
    end
    total_steps = nothing
    num_processors = nothing
    try
        total_steps = parse(Int, ARGS[2])
    catch e
        error("total_steps must be an integer.")
    end
    try
        num_processors = parse(Int, ARGS[3])
    catch e
        error("num_processors must be an integer.")
    end
    _run_parallel_marker_plotter(
        case_name      = case_name,
        total_steps    = total_steps,
        num_processors = num_processors
    )
    return nothing
end

"""
    _run_parallel_marker_plotter(
        case_name::String,
        total_steps::Int,
        num_processors::Int
    )::Nothing

Run the marker plotting script in parallel for a given range of time steps.

This function runs the Plot.jl script in parallel. Therefore, the Plot.jl
script must be in the same directory as the script that imports and executes
this function.

The plot command is of the form:

```bash
julia Plot.jl marker_plots case_name=<case_name> istart=<istart> iend=<iend>
```

This assumes that the command line plotter function `run_cl_plotter()` is being
passed the root path to the model output directory when called from the `Plot.jl` 
script.

# Arguments
- `case_name::String`: The name of the case.
- `total_steps::Int`: The total number of time steps.
- `num_processors::Int`: The number of processors to use.
"""
function _run_parallel_marker_plotter(;
    case_name::String,
    total_steps::Int,
    num_processors::Int
)::Nothing
    input_list = make_input_list(
        case_name      = case_name,
        total_steps    = total_steps,
        num_processors = num_processors
    )
    asyncmap(run_plot_script_wrapper, input_list; ntasks=num_processors)
    return nothing
end

"""
    make_input_list(
        case_name::String,
        total_steps::Int,
        num_processors::Int
    )::Vector{Vector{Union{Tuple{Int, Int}, String}}}

Make a list of input tuples for the parallel plotting loop.
"""
function make_input_list(;
    case_name::String,
    total_steps::Int,
    num_processors::Int
)::Vector{Vector{Union{Tuple{Int, Int}, String}}}
    steps_per_processor = total_steps ÷ num_processors
    ranges = make_time_step_ranges_list_for_each_processor(
        steps_per_processor=steps_per_processor,
        num_processors=num_processors,
        total_steps=total_steps
    )
    println(">> Plot parallel loop: steps_per_processor: $(steps_per_processor)")
    println(">> ranges: $(ranges)")
    input_list = Vector{Vector{Union{Tuple{Int, Int}, String}}}()
    for range_current in ranges
        push!(input_list, [range_current, case_name])
    end
    return input_list
end

"""
    make_time_step_ranges_list_for_each_processor(
        steps_per_processor::Int,
        num_processors::Int,
        total_steps::Int
    )::Vector{Tuple{Int, Int}}

Make a list of time step ranges for each processor.
"""
function make_time_step_ranges_list_for_each_processor(;
    steps_per_processor::Int,
    num_processors::Int,
    total_steps::Int
)::Vector{Tuple{Int, Int}}
    ranges = Vector{Tuple{Int, Int}}()
    for i in 0:(num_processors-1)
        push!(ranges, (i * steps_per_processor, (i + 1) * steps_per_processor - 1))
    end 
    if total_steps % num_processors != 0
        ranges[end] = (ranges[end][1], total_steps - 1)
    end
    return ranges
end

"""
    run_plot_script_wrapper(
        inputs::Vector{Union{Tuple{Int, Int}, String}}
    )::Nothing

Wrapper for running the plotting script.
"""
function run_plot_script_wrapper(inputs::Vector{Union{Tuple{Int, Int}, String}})::Nothing
    istart = inputs[1][1]
    iend = inputs[1][2]
    case_name = inputs[2]
    run_marker_plot_script(istart, iend, case_name)
    return nothing
end

"""
    run_marker_plot_script(
        istart::Int,
        iend::Int,
        case_name::String
    )::Nothing

Runs the plotting script for a given range of time steps using subprocess.
"""
function run_marker_plot_script(
    istart::Int,
    iend::Int,
    case_name::String
)::Nothing
    cmd = `julia Plot.jl marker_plots case_name=$(case_name) istart=$(istart) iend=$(iend)`
    try
        run(pipeline(cmd, stdout=devnull, stderr=devnull))
    catch e
        nothing
    end
    return nothing
end

end # module
