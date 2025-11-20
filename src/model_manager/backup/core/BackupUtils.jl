module BackupUtils

import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr

""" Check if object should be backed up.
"""
function check_for_backup(eb_obj::Any)::Bool
    check = false
    if !(eb_obj isa ParameterStr)
        if check_for_name(eb_obj)
            check = true
        end
    end
    return check
end

""" Check if object has attribute name.
"""
function check_for_name(eb_obj::Any)::Bool
    has_name = true
    try
        eb_obj.name
    catch
        has_name = false
    end
    return has_name
end

""" Make array for jld backup.
"""
function get_array(eb_obj)
    if eb_obj isa ParameterFloat
        ebarray = [eb_obj.value]
    elseif eb_obj isa ParameterInt
        ebarray = [eb_obj.value]
    else
        ebarray = eb_obj.array
    end
    return ebarray
end

end # module 