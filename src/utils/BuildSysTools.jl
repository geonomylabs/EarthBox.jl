module BuildSysTools

"""
SystemVectors is a struct that contains the vectors of the discretized system of 
equations.

# Attributes
- `Lii::Vector{Int64}`: 
    - Y-direction indices of basic grid for non-zero elements of the large
        matrix.

- `Ljj::Vector{Int64}`: 
    - X-direction indices of basic grid for non-zero matrix elements.

- `Li::Vector{Int64}`: 
    - Discretized large matrix row indices for non-zero matrix elements.

- `Lj::Vector{Int64}`: 
    - Discretized large matrix column indices for non-zero matrix elements.

- `Lv::Vector{Float64}`: 
    - Non-zero matrix values of the discretized system of equations.

"""
struct SystemVectors
    Lii::Vector{Int64}
    Ljj::Vector{Int64}
    Li::Vector{Int64}
    Lj::Vector{Int64}
    Lv::Vector{Float64}
    function SystemVectors(number_of_non_zero_elements::Int64)
        return new(
            zeros(Int64, number_of_non_zero_elements),
            zeros(Int64, number_of_non_zero_elements),
            zeros(Int64, number_of_non_zero_elements),
            zeros(Int64, number_of_non_zero_elements),
            zeros(Float64, number_of_non_zero_elements)
        )
    end
end

""" Clean non-zero matrix arrays.

# Arguments
- `nnz::Int64`: Number of non-zero matrix elements
- `Lii_tmp::Vector{Int64}`: Temporary basic grid row index for each non-zero matrix element 
- `Ljj_tmp::Vector{Int64}`: Temporary basic grid column index for each non-zero matrix element
- `Li_tmp::Vector{Int64}`: Temporary row index for each non-zero matrix element
- `Lj_tmp::Vector{Int64}`: Temporary column index for each non-zero matrix element
- `Lv_tmp::Vector{Float64}`: Temporary Non-zero matrix values

# Updated Arrays
- `Lii::Vector{Int64}`: Final basic grid row index for each non-zero matrix element 
- `Ljj::Vector{Int64}`: Final basic grid column index for each non-zero matrix element
- `Li::Vector{Int64}`: Final row index for each non-zero matrix element
- `Lj::Vector{Int64}`: Final column index for each non-zero matrix element
- `Lv::Vector{Float64}`: Final non-zero matrix values
"""
function clean_non_zero_arrays(
    nnz::Int64,
    Lii_tmp::Vector{Int64},
    Ljj_tmp::Vector{Int64},
    Li_tmp::Vector{Int64},
    Lj_tmp::Vector{Int64},
    Lv_tmp::Vector{Float64},
    Lii::Vector{Int64},
    Ljj::Vector{Int64},
    Li::Vector{Int64},
    Lj::Vector{Int64},
    Lv::Vector{Float64}
)
    for i in 1:nnz
        Lii[i] = Lii_tmp[i]
        Ljj[i] = Ljj_tmp[i]
        Li[i] = Li_tmp[i]
        Lj[i] = Lj_tmp[i]
        Lv[i] = Lv_tmp[i]
    end
end

""" Clean non-zero matrix arrays.

# Arguments
- `nnz::Int64`: Number of non-zero matrix elements
- `system_vectors_tmp::SystemVectors`: Temporary system vectors
- `system_vectors::SystemVectors`: Final system vectors

# Updated Arrays
- `system_vectors.Lii::Vector{Int64}`:
    - Final basic grid row index for each non-zero matrix element 
- `system_vectors.Ljj::Vector{Int64}`:
    - Final basic grid column index for each non-zero matrix element
- `system_vectors.Li::Vector{Int64}`:
    - Final row index for each non-zero matrix element
- `system_vectors.Lj::Vector{Int64}`:
    - Final column index for each non-zero matrix element
- `system_vectors.Lv::Vector{Float64}`:
    - Final non-zero matrix values

# Returns
- `system_vectors::SystemVectors`: Final system vectors
"""
function clean_non_zero_arrays!(
    nnz::Int64,
    system_vectors_tmp::SystemVectors
)::SystemVectors
    system_vectors = SystemVectors(nnz)
    for i in 1:nnz
        system_vectors.Lii[i] = system_vectors_tmp.Lii[i]
        system_vectors.Ljj[i] = system_vectors_tmp.Ljj[i]
        system_vectors.Li[i] = system_vectors_tmp.Li[i]
        system_vectors.Lj[i] = system_vectors_tmp.Lj[i]
        system_vectors.Lv[i] = system_vectors_tmp.Lv[i]
    end
    return system_vectors
end

end # module 