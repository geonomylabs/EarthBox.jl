""" Run the plotting script in parallel for a given range of time steps.

Usage:
    python plot_parallel.py case_name total_steps num_processors

    case_name : str
        The name of the case.

    total_steps : int
        The total number of time steps.

    num_processors : int
        The number of processors to use.
        
"""
import sys
from earthbox.plot_tools.utils.parallel import run_parallel_marker_plotter


if __name__ == "__main__":
    run_parallel_marker_plotter(
        case_name=sys.argv[1],
        total_steps=int(sys.argv[2]),
        num_processors=int(sys.argv[3])
        )
