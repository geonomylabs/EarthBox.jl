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

function main()
    run_earthbox(
        model_dir = "$(abspath(@__DIR__))",
        root_path_output_from_script = ROOT_PATH_OUTPUT
        )
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Runit.main()
end