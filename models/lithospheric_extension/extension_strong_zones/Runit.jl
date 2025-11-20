""" 
    Runit.jl

Run the EarthBox model and plot the results sequentially via command line 
execution in the background using the `run_earthbox` function. This script is 
required for remote execution. See documentation for `run_earthbox` for more 
details and usage examples.

!!! WARNING !!! This script is for Linux/Unix systems only.

# Quick Start

```bash
julia Runit.jl case_name=<case_name> run_model
```

```bash
julia Runit.jl case_name=<case_name> plot_scalars istart=1 iend=10
```

```bash
julia Runit.jl case_name=<case_name> plot_markers istart=1 iend=10
```

```bash
julia Runit.jl case_name=<case_name> run_model plot_markers plot_scalars istart=1 iend=10
```

"""
module Runit

using EarthBox
include("Model.jl")
import .Model: ROOT_PATH_OUTPUT

# Path to the EarthBox project used when running the model and plotting scripts.
const EB_PROJECT_PATH = "/home/$(get_username())/apps/julia/EarthBox"

function main()
    run_earthbox(
        model_dir = "$(abspath(@__DIR__))",
        # If the EarthBox project path is not provided, then it is assumed that 
        # EarthBox is installed in a standard Julia location or the JULIA_PROJECT 
        # environment variable is set equal to the EarthBox project path.
        # eb_project_path_from_script = EB_PROJECT_PATH,
        root_path_output_from_script = ROOT_PATH_OUTPUT
        )
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Runit.main()
end