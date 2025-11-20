module XdmfUtils

using Printf
import EarthBox.Arrays: ArrayUtils
import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D, AbstractEarthBoxArray2D

""" Convert integer to string with leading zeros.
"""
function intstr(number::Int)
    return @sprintf("%05d", number)
end

""" Get the xdmf number type of the array.
"""
function get_xdmf_number_type_for_array(
    array::Union{Array{Float64}, Array{Float32}, Array{Int64}, Array{Int16}}
)::String
    dtype_name = string(eltype(array))
    if occursin("Float", dtype_name)
        number_type = "Float"
    elseif occursin("Int", dtype_name)
        number_type = "Int"
    else
        error("Unable to determine the number type of array")
    end
    return number_type
end

function getoutform(
    obj::Union{AbstractEarthBoxArray1D, AbstractEarthBoxArray2D}
)::Union{Array{Float64}, Array{Int64}, Array{Int16}}
    return ArrayUtils.getoutform(obj)
end

end # module 