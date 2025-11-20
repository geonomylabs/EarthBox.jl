""" Module used to plot the characteristic new crust height vs model time

    Characteristic magmatic crust height is the height of the new crust formed
    at the mid-ocean ridge. This module extraction this information from an
    output log file called model_case#.log where # is the case number.

"""
import os
import re
import matplotlib.pyplot as plt
from earthbox.utilities.grep_tools import grep_file


def plot_magmatic_crust_height(
        output_directory_path: str,
        case_name: str,
        **kwargs
) -> None:
    """ Extract and plot data

    Inputs
    ----------
    output_directory_path : str
        Path to the directory containing log file containing data to be plotted.
        The log file should be named model_case#.log where # is the case number.

    case_name : str
        The case name of model run (e.g. case1).

    **kwargs
    ----------
    xdim : tuple
        Tuple containing the lower and upper limits of the x-axis.

    ydim : tuple
        Tuple containing the lower and upper limits of the y-axis.

    yticks : tuple
        Tuple containing the lower limit, upper limit, and interval of the y-axis
        ticks.

    """
    times, heights = extract_data(output_directory_path, case_name)
    plot_data(output_directory_path, case_name, times, heights, **kwargs)


def extract_data(
        output_directory_path: str,
        case_name: str
) -> tuple[list[float], list[float]]:
    """ Extract the characteristic new crust height and model time from log file
    
    re.findall
    ----------
    Searches data for all occurrences of the pattern using Python's re module
    and returns a list of tuples, where each tuple contains the captured 
    groups (crust height and model time).

    """
    file_path = f"{output_directory_path}/model_{case_name}.log"
    data = grep_file(file_path, ">> Characteristic new crust height")
    pattern = (
        r"Characteristic new crust height \(meters\): ([\d\.]+) model time \(Myr\): ([\d\.]+)"
        )
    matches = re.findall(pattern, data)
    heights = [float(match[0]) for match in matches]
    times = [float(match[1]) for match in matches]
    return times, heights


def plot_data(
        output_directory_path: str,
        case_name: str,
        times: list[float],
        heights: list[float],
        **kwargs
) -> None:
    """ Plot the characteristic new crust height vs model time
    """
    xdim = kwargs.get('xdim', (0, 10))
    ydim = kwargs.get('ydim', (0, 10000))
    yticks = kwargs.get('yticks', (0, 10000, 500))
    plt.figure(figsize=(10, 6))
    plt.plot(times, heights, linestyle='-', color='b')
    plt.title("Characteristic New Crust Height vs Model Time")
    plt.xlabel("Model Time (Myr)")
    plt.ylabel("Characteristic New Crust Height (meters)")
    plt.xlim(xdim[0], xdim[1])
    plt.ylim(ydim[0], ydim[1])
    plt.yticks(range(yticks[0], yticks[1], yticks[2]))
    plt.grid(True)
    plots_directory_path = os.path.join(output_directory_path, 'plots')
    if not os.path.exists(plots_directory_path):
        os.makedirs(plots_directory_path)
    file_path = os.path.join(
        plots_directory_path,
        f"magmatic_crust_height_{case_name}.png"
        )
    plt.savefig(file_path)
