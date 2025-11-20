""" Execute a script on a remote machine in the background.
"""
from earthbox.utils.run_tools import model_loop


def run_models():
    """ Run models from model groups in background on remote hosts.

    This script runs the runit.py scripts stored in the model directories on 
    remote hosts. This script assumes that host names and ssh keys are set up 
    properly and that the paths to output directories defined the the path.py 
    script in model directories have paths that point to existing directories on
    nodes.

    Parameters
    ----------
    models : dict[str, list[str]]
        A dictionary of model info and corresponding host names. The models 
        dictionary should have the following format:

        models = {
            "node1": {
                "model_base_name": "model_base_name1",
                "model_case_names": ["case0", "case1"]
                },
            "node2": {
                "model_base_name": "model_base_name2",
                "model_case_names": ["case0", "case1"]
                }
            }

        where "node1" and "node2 are the host names, "model_base_name" is the 
        model directory name that contains the model scripts, and 
        "model_case_names" is a list of model case names as defined in the
        model scripts. 
         
    model_group_dir_path : str
        The path to the directory that contains the model directories. Model
        directories should have names equal to the value of the key 
        "model_base_name".

    In order for this to work, comment out the following from .bashrc if it exists:
        case $- in
            *i*) ;;
            *) return;;
        esac

    """
    model_group_name = 'extension_to_sfs'

    base_names = {
        "extension_to_sfs": ["ssatl"],
        }

    model_group_paths = {
        "extension_to_sfs": r"/home/ekneller/apps/python/earthbox/examples/extension_to_sfs",
        }

    models = {
        #"peridotite": {
        #    'model_base_name': base_names[model_group_name][0],
        #    'model_case_names': ['case0', 'case1', 'case2', 'case3']
        #    },
        "eclogite": {
            'model_base_name': base_names[model_group_name][0],
            'model_case_names': ['case4', 'case5', 'case6', 'case7']
            },
        #"lherzolite": {
        #    'model_base_name': base_names[model_group_name][0],
        #    'model_case_names': ['case4', 'case5', 'case6', 'case7']
        #    },
        #"plagioclase": {
        #    'model_base_name': base_names[model_group_name][0],
        #    'model_case_names': ['case8', 'case9', 'case10', 'case11']
        #    },
        #"chromite": {
        #    'model_base_name': base_names[model_group_name][0],
        #    'model_case_names': ['case12', 'case13', 'case14', 'case15']
        #    },
        #"pyroxene": {
        #    'model_base_name': base_names[model_group_name][0],
        #    'model_case_names': ['case16', 'case17', 'case18', 'case19']
        #    },
        #"dunite": {
        #    'model_base_name': base_names[model_group_name][0],
        #    'model_case_names': ['case20', 'case21', 'case22', 'case23']
        #    },
        #"komatiite": {
        #    'model_base_name': base_names[model_group_name][0],
        #    'model_case_names': ['case24']
        #    },
        #"spinel": {
        #     'model_base_name': base_names[model_group_name][0],
        #     'model_case_names': ['case42', 'case43', 'case44', 'case47']
        #     }
        }

    model_loop(models, model_group_paths[model_group_name])


if __name__ == "__main__":
    run_models()
