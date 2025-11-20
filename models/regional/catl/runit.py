""" Run the model and copy the output to the output directory.
"""
from paths import get_output_path
from earthbox.utils.run_tools import execute_earthbox_script


def run_earthbox() -> None:
    """ Run the model from the model directory.
    """
    execute_earthbox_script(
        command_type="model",
        model_logfile_name="model.log",
        final_log_dir_path=get_output_path()
        )
    execute_earthbox_script(command_type="plot_markers")


if __name__ == "__main__":
    run_earthbox()
