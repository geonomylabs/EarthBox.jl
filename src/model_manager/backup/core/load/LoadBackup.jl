module LoadBackup

include("SetObjects.jl")

using Printf
import JLD2
import EarthBox.EarthBoxDtypes: ObjDictType
import .SetObjects: setobj
import ..BackupUtils: check_for_backup

"""
    load_backup_jld2(obj_dict::Dict{String,Any}, output_dir::String)::Nothing

Load arrays from model backup file and set earthbox objects.

# Arguments
- `obj_dict::ObjDictType`: Dictionary of earthbox parameter and array objects to update.
  These objects have a standard set of attributes like name, description, and value for scalar
  parameters or array for arrays. Note that both parameters and arrays are stored in the jld
  file as arrays. To get the value of a scalar from an jld array just use the first element.
- `output_dir::String`: Directory containing model backup file.

# Throws
- `SystemError`: If the model backup file is not found
"""
function load_backup_jld2(obj_dict::ObjDictType, output_dir::String)::Nothing
    jld_filename = joinpath(output_dir, "model_backup.jld")
    if isfile(jld_filename)
        JLD2.jldopen(jld_filename, "r") do file
            for (eb_obj_name, eb_obj) in obj_dict
                if check_for_backup(eb_obj)
                    if haskey(file, eb_obj_name)
                        data = file[eb_obj_name]
                        setobj(eb_obj, data)
                    else
                        @printf(">> No data found in jld file for %s\n", eb_obj_name)
                    end
                end
            end
        end
    else
        error(">> No model backup file found in $output_dir. " *
              "Set restart_from_backup to false when initiating Earthbox to run model.")
    end
    
    return nothing
end

end # module 