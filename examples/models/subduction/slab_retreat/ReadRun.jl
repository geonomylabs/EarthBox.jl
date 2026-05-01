"""
   ReadRun.jl

Run EarthBox model using yaml input files and paths defined by functions.

Usage from command line:

```bash
% julia ReadRun.jl
```

Usage from REPL:

```julia
include("ReadRun.jl")
ReadRun.read_run()
```
"""
module ReadRun

using EarthBox

const MODEL_PATH = joinpath(@__DIR__, "model.yml")
const MATERIAL_MODEL_PATH = joinpath(@__DIR__, "materials.yml")
const MATERIAL_COLLECTION = MaterialLibrary().gerya2019
const MODEL_OUTPUT_PATH = "/mnt/extradrive1/earthbox_output/slab_retreat_output"

function read_run()::Nothing
    marker_output_from_user = Dict{String, Bool}(
        "marker_age" => true,
        "marker_rho" => true,
        "marker_serpentinization" => false,
        "marker_extractable_meltfrac" => false,
        "marker_extracted_meltfrac" => false,
        "marker_strain_rate_plastic" => true,
    )
    run_model(
        make_backup             = false, 
        restart_from_backup     = false, 
        use_mumps               = false, 
        use_mumps_internal      = true, 
        nprocs                  = 1, 
        marker_output_from_user = marker_output_from_user,
        model_input_file        = MODEL_PATH, 
        materials_input_file    = MATERIAL_MODEL_PATH, 
        materials_library_file  = MATERIAL_COLLECTION.path, 
        model_output_path       = MODEL_OUTPUT_PATH
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