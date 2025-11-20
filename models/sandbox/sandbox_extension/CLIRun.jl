"""
    CLRun.jl

Run EarthBox model using yamal input files and command line arguments.

Usage:

```bash
% julia CLIRun.jl model_input_file="model.yml" materials_input_file="materials.yml" materials_library_file="/path/to/materials_collection_file.yml" model_output_path="path/to/model_output marker_output='{marker_x=false, marker_y=true}' 
```

"""
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
