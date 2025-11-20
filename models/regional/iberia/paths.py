""" Define paths to material library and output directory.
"""
import os


def get_output_path():
    """ Define the path to the output directory.
    """
    drive_number = 1
    base_path = r'/mnt/extradrive' + str(drive_number) + r'/earthbox_output/extension_to_sfs'
    model_name = 'iberia_case0'
    output_dir_name = model_name + '_output'
    if not os.path.exists(base_path):
        os.makedirs(base_path)
    return os.path.join(base_path, output_dir_name)
