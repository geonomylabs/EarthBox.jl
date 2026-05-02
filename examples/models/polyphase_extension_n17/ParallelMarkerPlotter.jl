""" 
    ParallelMarkerPlotter.jl

Run the marker plotting script `Plot.jl` in parallel.

Usage:
    julia ParallelMarkerPlotter.jl <case_name> <total_steps> <num_processors>

where:
- <case_name>: The name of the case with the format case# where # is an integer.
- <total_steps>: The total number of time steps.
- <num_processors>: The number of processors to use.
        
"""
module ParallelMarkerPlotter

using EarthBox

function main()::Nothing
    run_parallel_marker_plotter()
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ParallelMarkerPlotter.main()
end

