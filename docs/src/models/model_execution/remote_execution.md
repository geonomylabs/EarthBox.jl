# Batch Remote Execution

`Runit.jl` scripts that call the function [EarthBox.run_earthbox](@ref) to run 
`Model.jl` and `Plot.jl scripts` in the background via command line can be run 
on remote machines using the function [`remote_model_loop`](@ref "Remote Model Loop"). See
[`execute_remote_script_in_background`](@ref "Executing Remote Scripts in the Background") for details how passwordless `ssh` accessed is used for remote execution. 

!!! warning "Linux/Unix Only"
    Remote model execution functions are compatible with Linux/Unix operating 
    systems only.

```julia
module RunitRemote

using EarthBox

function run_remote_models()
    model_dir_path = "/mnt/extradrive1/nfs/apps/julia/EarthBox/models/extension_to_sfs/extension_strong_zones"
    models = Dict(
        "hostname1" => ["case0", "case1"],
        "hostname2" => ["case0", "case1", "case2"],
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
```
## Remote Model Loop

```@docs
EarthBox.RunTools.remote_model_loop
```

## Executing Remote Scripts in the Background

```@docs
EarthBox.RunTools.execute_remote_script_in_background
```