""" Run the model and copy the output to the output directory.
"""
from paths import get_output_path
from earthbox.utils.get_paths import get_model_case_name
from earthbox.utils.run_tools import run_earthbox


if __name__ == "__main__":
    model_case_name = get_model_case_name()
    run_earthbox(
        model_case_name=model_case_name,
        final_log_dir_path=get_output_path(
            model_case_name=model_case_name, drive_number_id=2),
        run_model=True,
        plot_markers=True,
        plot_scalars=False,
        plot_magmatic_crust_height=False
        )
