module BuildSysTools

"""
SystemVectors holds the preallocated build buffers and downstream-facing
output buffers for the discretized system of equations.

# Build buffers (sized to `nnz_max`, the worst-case non-zero count)
- `Lii::Vector{Int64}`: Y-direction (basic-grid row) index of each
    non-zero matrix element. Diagnostic; not consumed downstream.
- `Ljj::Vector{Int64}`: X-direction (basic-grid column) index of each
    non-zero matrix element. Diagnostic; not consumed downstream.
- `Li::Vector{Int64}`: Discretized large-matrix row index for each
    non-zero matrix element. Build pipeline writes here using an `inz`
    cursor.
- `Lj::Vector{Int64}`: Discretized large-matrix column index. Same.
- `Lv::Vector{Float64}`: Non-zero matrix values. Same.

# Output buffers (resized to `nnz` by `clean_non_zero_arrays!`)
- `Li_out::Vector{Int64}`: Length-`nnz` packed prefix of `Li` after each
    build. Capacity preserved at `nnz_max` so `resize!`/`copyto!` cycles
    do not allocate.
- `Lj_out::Vector{Int64}`: Same for `Lj`.
- `Lv_out::Vector{Float64}`: Same for `Lv`.

Downstream consumers (`sparse(...)`, MUMPS) read only the `*_out` fields,
which behave as concrete `Vector` of length `nnz`. `Lii`/`Ljj` are not
mirrored because no downstream consumer reads them.
"""
struct SystemVectors
    Lii::Vector{Int64}
    Ljj::Vector{Int64}
    Li::Vector{Int64}
    Lj::Vector{Int64}
    Lv::Vector{Float64}
    Li_out::Vector{Int64}
    Lj_out::Vector{Int64}
    Lv_out::Vector{Float64}
    function SystemVectors(number_of_non_zero_elements::Int64)
        n = number_of_non_zero_elements
        return new(
            zeros(Int64, n),
            zeros(Int64, n),
            zeros(Int64, n),
            zeros(Int64, n),
            zeros(Float64, n),
            zeros(Int64, n),
            zeros(Int64, n),
            zeros(Float64, n)
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
    sv::SystemVectors
)::SystemVectors
    # Mirror the first `nnz` entries of the build buffers into the
    # downstream-facing output buffers without allocating: resize! to nnz
    # (capacity preserved at nnz_max), then copyto! the prefix.
    resize!(sv.Li_out, nnz)
    resize!(sv.Lj_out, nnz)
    resize!(sv.Lv_out, nnz)
    copyto!(sv.Li_out, 1, sv.Li, 1, nnz)
    copyto!(sv.Lj_out, 1, sv.Lj, 1, nnz)
    copyto!(sv.Lv_out, 1, sv.Lv, 1, nnz)
    return sv
end

end # module 