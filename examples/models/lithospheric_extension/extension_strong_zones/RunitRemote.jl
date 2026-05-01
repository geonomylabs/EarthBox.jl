module RunitRemote

using EarthBox

"""
    run_remote_models()

Run models from model groups in background on remote hosts.

!!! WARNING !!! This script is for Linux/Unix systems only.

This script runs the Runit.jl script stored in the model directories on 
remote hosts. This script assumes that host names and ssh keys are set up 
properly and that the paths to output directories defined the the path.py 
script in model directories have paths that point to existing directories on
nodes.

In order for this to work, comment out the following from .bashrc if it exists:
```bash
case \$- in
    *i*) ;;
    *) return;;
esac
```
"""
function run_remote_models()
    # Path to the model directory where the Runit.jl script is stored.
    # Model.jl and Plot.jl should also be in this directory.
    model_dir_path = "/mnt/extradrive1/nfs/apps/julia/EarthBox/models/extension_to_sfs/extension_strong_zones"
    # Dictionary of host names and model case names that will be run on the 
    # remote hosts.
    models = Dict(
        "peridotite" => ["case0", "case1"],
        #"eclogite" => ["case0", "case1", "case2", "case3"],
        #"lherzolite" => ["case4", "case5", "case6", "case7"],
        #"plagioclase" => ["case8", "case9", "case10", "case11"],
        #"chromite" => ["case12", "case13", "case14", "case15"],
        #"pyroxene" => ["case18"],
        #"dunite" => ["case20", "case21", "case22", "case23"],
        #"komatiite" => ["case24"],
        #"spinel" => ["case42", "case43", "case44", "case47"]
    )
    remote_model_loop(
        models,
        model_dir_path,
        run_model = true,
        plot_markers = false,
        plot_scalars = false,
        istart = 1,
        iend = 100
    )
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    RunitRemote.run_remote_models()
end
