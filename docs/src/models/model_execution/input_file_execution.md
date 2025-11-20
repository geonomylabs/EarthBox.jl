# Executing Models Using Input Files

The `run_model` function can be used in a script with arguments based on user 
defined functions or via command line.

```@docs
EarthBox.run_model
```

## Execution Using `run_model` and Command-line Arguments (CLIRun.jl)

A command-line based script, typically named `CLIRun.jl`, can be developed using 
the `run_model` function and functions from the `GetArgs`:

```@julia
module CLIRun

using EarthBox

function cl_run()
    run_model(
       make_backup             = GetArgs.get_make_backup_from_args(), 
       restart_from_backup     = GetArgs.get_restart_from_backup_from_args(), 
       use_mumps               = GetArgs.get_use_mumps_from_args(), 
       use_mumps_internal      = GetArgs.get_use_mumps_internal_from_args(), 
       nprocs                  = GetArgs.get_nprocs_from_args(), 
       marker_output_from_user = GetArgs.get_marker_output_dict_from_args(),
       model_input_file        = GetArgs.get_model_input_file_from_args(), 
       materials_input_file    = GetArgs.get_materials_input_file_from_args(), 
       materials_library_file  = GetArgs.get_materials_library_file_from_args(), 
       model_output_path       = GetArgs.get_model_output_path_from_args()
    )
    return nothing
end

function main()::Nothing
    cl_run()
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    CLIRun.main()
end
```

If `model.yml` and `materials.yml` are in the same directory as `CLRun.jl` then
the model can be run with the following command:

```
% julia CLIRun.jl model_input_file="model.yml" materials_input_file="materials.yml" materials_library_file="/path/to/materials_collection_file.yml" model_output_path="path/to/model_output marker_output='{marker_x=false, marker_y=true}'" 
```

!!! tip "`model_output_path`"
    If you omit the `model_output_path` argument output will be sent to a directory 
    called `output` in the current working directory.

!!! tip "material collection files"
    Material collection files are stored in src/material_library/registries.

!!! warning "marker_output"
    When defining the marker output flag dictionary via command line ensure that
    either no spaces exist on either side of the = sign or use a space on both 
    sides of the = sign. For example, the following are valid command-line 
    arguments:
    - `marker_output={marker_x=true, marker_y=true}`
    - `marker_output = {marker_x=true, marker_y=true}`


## Execution Using `run_model` and User-defined Functions (ReadRun.jl) 

```julia
module ReadRun

using EarthBox

const MODEL_PATH = joinpath(@__DIR__, "model.yml")
const MATERIAL_MODEL_PATH = joinpath(@__DIR__, "materials.yml")
const MATERIAL_COLLECTION = MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1
const MODEL_OUTPUT_PATH = "/mnt/extradrive1/earthbox_output/extension_to_sfs/sdr_formation_output"

function read_run()::Nothing
    marker_output_from_user = Dict{String, Bool}(
        "marker_age" => true,
        "marker_rho" => true,
        "marker_serpentinization" => true,
        "marker_extractable_meltfrac" => false,
        "marker_extracted_meltfrac" => false,
        "marker_strain_rate_plastic" => true,
    )
    run_model(
        make_backup            = false, 
        restart_from_backup    = false, 
        use_mumps              = false, 
        use_mumps_internal     = true, 
        nprocs                 = 1, 
        marker_output_from_user = marker_output_from_user,
        model_input_file       = MODEL_PATH, 
        materials_input_file   = MATERIAL_MODEL_PATH, 
        materials_library_file = MATERIAL_COLLECTION.path, 
        model_output_path      = MODEL_OUTPUT_PATH
        )
    return nothing
end

function main()::Nothing
    read_run()
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ReadRun.main()
end
```

The script above can be executed via command line:

```bash
%julia ReadRun.jl
```

or from the REPL,

```julia
include("ReadRun.jl")
ReadRun.read_un()
```