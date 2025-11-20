module ArrayUtils

import LinearAlgebra
import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D, AbstractEarthBoxArray2D

mutable struct OutputFormat
    fac1::Float64
    fac2::Float64
    units::String
    name::String
    log10::Bool
    fname::String
    header::String
end

function setzeros!(eb_array::AbstractEarthBoxArray1D)
    eb_array.array .= 0.0
    return eb_array.array
end

function setzeros!(eb_array::AbstractEarthBoxArray2D)
    fill!(eb_array.array, 0.0)
    return eb_array.array
end

function print_min_max(name::String, eb_array::Vector{Float64})
    min_val = LinearAlgebra.minimum(eb_array)
    max_val = LinearAlgebra.maximum(eb_array)
    mean_val = sum(eb_array) / length(eb_array)
    println("$(name) : min : $min_val : max : $max_val : mean : $mean_val")
end

function print_min_max(name::String, eb_array::Matrix{Float64})
    min_val = LinearAlgebra.minimum(eb_array)
    max_val = LinearAlgebra.maximum(eb_array)
    mean_val = sum(eb_array) / length(eb_array)
    println("$(name) : min : $min_val : max : $max_val : mean : $mean_val")
end

function print_min_max(name::String, eb_array::Vector{Int64})
    min_val = LinearAlgebra.minimum(eb_array)
    max_val = LinearAlgebra.maximum(eb_array)
    mean_val = sum(eb_array) / length(eb_array)
    println("$(name) : min : $min_val : max : $max_val : mean : $mean_val")
end

function print_min_max(eb_array::AbstractEarthBoxArray1D)
    min_val = LinearAlgebra.minimum(eb_array.array)
    max_val = LinearAlgebra.maximum(eb_array.array)
    mean_val = sum(eb_array.array) / length(eb_array.array)
    println("$(eb_array.name) : min : $min_val : max : $max_val : mean : $mean_val")
end

function print_min_max(eb_array::AbstractEarthBoxArray2D)
    min_val = LinearAlgebra.minimum(eb_array.array)
    max_val = LinearAlgebra.maximum(eb_array.array)
    mean_val = sum(eb_array.array) / length(eb_array.array)
    println("$(eb_array.name) : min : $min_val : max : $max_val : mean : $mean_val")
end

function getoutform(eb_array::AbstractEarthBoxArray1D)::Union{Vector{Float64}, Vector{Int64}, Vector{Int16}}
    fac1 = eb_array.outform.fac1
    fac2 = eb_array.outform.fac2
    use_log10 = eb_array.outform.log10

    data = eb_array.array .* fac1 .+ fac2
    if use_log10
        return log10.(abs.(data))
    else
        return data
    end
end

function getoutform(eb_array::AbstractEarthBoxArray2D)::Union{Matrix{Float64}, Matrix{Int64}, Matrix{Int16}}
    fac1 = eb_array.outform.fac1
    fac2 = eb_array.outform.fac2
    use_log10 = eb_array.outform.log10

    data = eb_array.array .* fac1 .+ fac2
    if use_log10
        return log10.(abs.(data))
    else
        return data
    end
end

function check_limits(
    msg::String,
    eb_array::Union{Vector{Float64}, Vector{Int64}, Matrix{Float64}, Matrix{Int64}}
)::Nothing
    println(">> $msg : $(LinearAlgebra.minimum(eb_array)) : $(LinearAlgebra.maximum(eb_array))")
end

end # module
