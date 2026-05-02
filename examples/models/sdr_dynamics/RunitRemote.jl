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

If `using EarthBox` fails on a node only for **remote** runs, non-interactive SSH
often does not see `JULIA_DEPOT_PATH` from `.bashrc` (exports after the early
return above are skipped). Move `JULIA_DEPOT_PATH` exports to the **top** of 
`~/.bashrc` before the `case \$-` block.
"""
function run_remote_models()
    # Path to the model directory where the Runit.jl script is stored.
    # Model.jl and Plot.jl should also be in this directory.
    model_dir_path = "/home/$(get_username())/apps/julia/EarthBox_models/sdr_dynamics"
    # Dictionary of host names and model case names that will be run on the 
    # remote hosts.
    models = Dict(
        #"peridotite" => ["case0", "case1"],
        #"eclogite" => ["case2", "case3"],
        "lherzolite" => ["case4", "case5"],
        #"plagioclase" => ["case6", "case7"],
        #"chromite" => ["case8", "case9"],
        #"pyroxene" => ["case10","case11"],
	    #"spinel" => ["case12", "case13", "case23"],
        #"komatiite" => ["case14", "case15"],
        #"dunite" => ["case16", "case17"],
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
