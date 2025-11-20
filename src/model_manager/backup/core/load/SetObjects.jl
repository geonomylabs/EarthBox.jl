module SetObjects

import EarthBox.Arrays: ArrayTypes
import EarthBox: Parameters

function setobj(eb_obj::Parameters.ParameterInt, array::Vector{Int64})::Nothing
    eb_obj.value = array[1]
    return nothing
end

function setobj(eb_obj::Parameters.ParameterFloat, array::Vector{Float64})::Nothing
    eb_obj.value = array[1]
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.InternalBcArrayFloat.InternalBcArrayFloatState,
    array::Vector{Float64}
)::Nothing
    check_1d_float_arrays(array, eb_obj.array)
    for i in 1:length(array)
        eb_obj.array[i] = array[i]
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.InternalBcArrayInt.InternalBcArrayIntState,
    array::Vector{Int64}
)::Nothing
    check_1d_int_arrays(array, eb_obj.array)
    for i in 1:length(array)
        eb_obj.array[i] = array[i]
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.BcArrayFloat.BcArrayFloatState,
    array::Matrix{Float64}
)::Nothing
    check_2d_float_arrays(array, eb_obj.array)
    for i in 1:size(array, 1)
        for j in 1:size(array, 2)
            eb_obj.array[i, j] = array[i, j]
        end
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.CarbArray2D.CarbArray2DState,
    array::Matrix{Float64}
)::Nothing
    check_2d_float_arrays(array, eb_obj.array)
    for i in 1:size(array, 1)
        for j in 1:size(array, 2)
            eb_obj.array[i, j] = array[i, j]
        end
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.TopoArray2D.TopoArray2DState,
    array::Matrix{Float64}
)::Nothing
    check_2d_float_arrays(array, eb_obj.array)
    for i in 1:size(array, 1)
        for j in 1:size(array, 2)
            eb_obj.array[i, j] = array[i, j]
        end
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.GridArray1D.GridArray1DState,
    array::Vector{Float64}
)::Nothing
    check_1d_float_arrays(array, eb_obj.array)
    for i in 1:length(array)
        eb_obj.array[i] = array[i]
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.ScalarArray2D.ScalarArray2DState,
    array::Matrix{Float64}
)::Nothing
    check_2d_float_arrays(array, eb_obj.array)
    for i in 1:size(array, 1)
        for j in 1:size(array, 2)
            eb_obj.array[i, j] = array[i, j]
        end
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.SolutionArray1D.SolutionArray1DState,
    array::Vector{Float64}
)::Nothing
    check_1d_float_arrays(array, eb_obj.array)
    for i in 1:length(array)
        eb_obj.array[i] = array[i]
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.RhsHeatArray1D.RhsHeatArray1DState,
    array::Vector{Float64}
)::Nothing
    check_1d_float_arrays(array, eb_obj.array)
    for i in 1:length(array)
        eb_obj.array[i] = array[i]
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.MarkerArrayFloat1D.MarkerArrayFloat1DState{T},
    array::Vector{T}
)::Nothing where {T<:AbstractFloat}
    check_1d_float_arrays(array, eb_obj.array)
    for i in 1:length(array)
        eb_obj.array[i] = array[i]
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.MarkerArrayInt1D.MarkerArrayInt1DState{T},
    array::Vector{T}
)::Nothing where {T<:Integer}
    check_1d_int_arrays(array, eb_obj.array)
    for i in 1:length(array)
        eb_obj.array[i] = array[i]
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.MaterialArrayFloat1D.MaterialArrayFloat1DState,
    array::Vector{Float64}
)::Nothing
    check_1d_float_arrays(array, eb_obj.array)
    for i in 1:length(array)
        eb_obj.array[i] = array[i]
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.MaterialArrayFloat2D.MaterialArrayFloat2DState,
    array::Matrix{Float64}
)::Nothing
    check_2d_float_arrays(array, eb_obj.array)
    for i in 1:size(array, 1)
        for j in 1:size(array, 2)
            eb_obj.array[i, j] = array[i, j]
        end
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.RhsStokesArray1D.RhsStokesArray1DState,
    array::Vector{Float64}
)::Nothing
    check_1d_float_arrays(array, eb_obj.array)
    for i in 1:length(array)
        eb_obj.array[i] = array[i]
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.MaterialArrayInt1D.MaterialArrayInt1DState,
    array::Vector{Int64}
)::Nothing
    check_1d_int_arrays(array, eb_obj.array)
    for i in 1:length(array)
        eb_obj.array[i] = array[i]
    end
    return nothing
end

function setobj(
    eb_obj::ArrayTypes.MaterialArrayInt2D.MaterialArrayInt2DState,
    array::Matrix{Int64}
)::Nothing
    check_2d_int_arrays(array, eb_obj.array)
    for i in 1:size(array, 1)
        for j in 1:size(array, 2)
            eb_obj.array[i, j] = array[i, j]
        end
    end
    return nothing
end

function check_1d_float_arrays(
    array_jld::Union{Array{Float64}, Array{Float32}},
    array_eb::Union{Array{Float64}, Array{Float32}}
)::Nothing
    if length(array_jld) != length(array_eb)
        throw(ErrorException("1D float array sizes do not match."))
    end
    return nothing
end

function check_2d_float_arrays(
    array_jld::Union{Array{Float64}, Array{Float32}},
    array_eb::Union{Array{Float64}, Array{Float32}}
)::Nothing
    if size(array_jld, 1) != size(array_eb, 1)
        throw(ErrorException("2D float array [1] sizes do not match."))
    end
    if size(array_jld, 2) != size(array_eb, 2)
        throw(ErrorException("2D float array [2] sizes do not match."))
    end
    return nothing
end

function check_1d_int_arrays(
    array_jld::Array{T},
    array_eb::Array{T}
)::Nothing where {T<:Integer}
    if length(array_jld) != length(array_eb)
        throw(ErrorException("1D int array sizes do not match."))
    end
    return nothing
end

function check_2d_int_arrays(
    array_jld::Array{Int64},
    array_eb::Array{Int64}
)::Nothing
    if size(array_jld, 1) != size(array_eb, 1)
        throw(ErrorException("2D int array [1] sizes do not match."))
    end
    if size(array_jld, 2) != size(array_eb, 2)
        throw(ErrorException("2D int array [2] sizes do not match."))
    end
    return nothing
end

end # module 