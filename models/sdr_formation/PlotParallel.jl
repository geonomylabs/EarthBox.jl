""" Run the plotting script in parallel for a given range of time steps.

Usage:
    julia PlotParallel.jl case_name total_steps num_processors

    case_name : String
        The name of the case.

    total_steps : Int
        The total number of time steps.

    num_processors : Int
        The number of processors to use.
        
"""
module PlotParallel

using EarthBox

export main

function main()::Nothing
    if length(ARGS) != 3
        error("Usage: julia PlotParallel.jl case_name total_steps num_processors")
    end
    
    case_name = ARGS[1]
    total_steps = parse(Int, ARGS[2])
    num_processors = parse(Int, ARGS[3])
    
    run_parallel_marker_plotter(
        case_name=case_name,
        total_steps=total_steps,
        num_processors=num_processors
    )
end

# Run main function if this script is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

end # module
