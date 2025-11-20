module MakeBackup

using Printf
import JLD2
import EarthBox.EarthBoxDtypes: ObjDictType
import ..BackupUtils: get_array, check_for_backup

"""
    make_backup_jld2(obj_dict::ObjDictType, output_dir::String)::Nothing

Make jld2 file for model backup.

# Arguments
- `obj_dict::ObjDictType`: Dictionary of earthbox parameter and array objects to backup.
  These objects have a standard set of attributes like name, description, and value for scalar
  parameters or array for arrays.
- `output_dir::String`: Directory to save model backup file.

# Throws
- `SystemError`: If unable to create or write to the backup file
"""
function make_backup_jld2(obj_dict::ObjDictType, output_dir::String)::Nothing
    jld_filename = joinpath(output_dir, "model_backup.jld")
    
    JLD2.jldopen(jld_filename, "w") do file
        for (eb_obj_name, eb_obj) in obj_dict
            if check_for_backup(eb_obj)
                ebarray = get_array(eb_obj)
                file[eb_obj_name] = ebarray
            end
        end
    end
    
    return nothing
end

end # module 