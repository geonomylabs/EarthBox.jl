""" Define paths to material library and output directory.
"""
from earthbox.utils.get_paths import get_output_path as _get_output_path


def get_output_path(**kwargs) -> str:
    """ Define the path to the output directory.

    Inputs
    ------
    kwargs: dict, optional
    - drive_root_path: str
        The root path to the drive such a /mnt if drive is mounted at /mnt/drive_name.

    - drive_base_name: str
        The base name of the drive where drives are named as follows:
            drive_base_name + drive_number

    - drive_number_id: int
        Integer ID of the drive.

    - local_group_path_on_drive: str
        The local path to output group directory on the drive where a group is
        a collection of models of similar type.

    - model_base_name: str
        The base name of the model.

    - model_case_name: str
        The case name of the model.

    Returns
    -------
    output_path: str
        The path to the output directory such as 
        /mnt/extradrive1/earthbox_output/extension_to_sfs/constant_rate/case0
        where :
            - drive_root_path = '/mnt'
            - drive_base_name = 'extradrive'
            - drive_number = 1
            - local_group_path_on_drive = 'earthbox_output/extension_to_sfs'
            - model_base_name = 'constant_rate'
            - model_case_name = 'case0'

    """
    return _get_output_path(
        drive_root_path=kwargs.get(
            'drive_root_path', r'/mnt'),
        drive_base_name=kwargs.get(
            'drive_base_name', 'extradrive'),
        drive_number_id=kwargs.get(
            'drive_number_id', 1),
        local_group_path_on_drive=kwargs.get(
            'local_group_path_on_drive', 'earthbox_output/extension_to_sfs'),
        model_base_name=kwargs.get(
            'model_base_name', 'ssatl'),
        model_case_name=kwargs.get(
            'model_case_name', 'case0')
        )
