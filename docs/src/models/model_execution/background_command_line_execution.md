# Background Command-line Execution

The function [EarthBox.run_earthbox](@ref) enables sequential model running and 
post-processing via background command-line execution on Linux/Unix systems. 

In order to have compatibility with remote model execution scripts the 
[EarthBox.run_earthbox](@ref) function should be used inside a script called 
`Runit.jl` that can be executed via command line and is located in the model 
directory with [`Model.jl`](@ref) and [`Plot.jl`](@ref) scripts. 

# `Runit.jl`.

Creating a script called `Runit.jl` in the model directory that calls the 
function [EarthBox.run_earthbox](@ref) allows the use to run models and 
post-processing functions sequentially. See [EarthBox.run_earthbox](@ref) for 
command-line arguments that are used by [EarthBox.run_earthbox](@ref). 

!!! note "remote and local batch execution"
    Scripts named `Runit.jl` in a model directory scripts are executed by 
    [`remote_model_loop`](@ref "Remote Model Loop") for remote execution and 
    by [`local_model_loop`](@ref "Local Model Loop") for local batch execution. 

!!! warning "Runit.jl file name"
    Since remote and local batch processing functions look for a script called 
    `Runit.jl` in the model directory it is critical to use the name `Runit.jl`.

!!! warning "local of Runit.jl scripts"
    `Runit.jl` scripts must be located in a model directory along with ['Model.jl']
    and ['Plot.jl'] scripts.

The following is an example of a `Runit.jl` script that would enable sequential
execution of models and post-processing functions in addition to remote and local
batch execution.

```julia
module Runit

using EarthBox
include("Model.jl")

# Get the model output root path constants from Model.jl
import .Model: ROOT_PATH_OUTPUT

# Path to the EarthBox project used when running the model and plotting scripts.
const EB_PROJECT_PATH = "/home/$(get_username())/apps/julia/EarthBox"

function main()
    run_earthbox(
        model_dir = "$(abspath(@__DIR__))",
        # If the EarthBox project path is not provided, then it is assumed that 
        # EarthBox is installed in a standard Julia location or the JULIA_PROJECT 
        # environment variable is set equal to the EarthBox project path.
        eb_project_path_from_script = EB_PROJECT_PATH,
        root_path_output_from_script = ROOT_PATH_OUTPUT
        )
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Runit.main()
end
```

```@docs
EarthBox.run_earthbox
```

```@docs
EarthBox.execute_earthbox_script
```


