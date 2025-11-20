module SystemReader

using JLD2

function read_system_of_equations_binary(filepath::String)
    open(filepath, "r") do io
        N = read(io, Int64)
        array_length = read(io, Int64)
        Li = Vector{Int64}(undef, array_length)
        read!(io, Li)
        Lj = Vector{Int64}(undef, array_length)
        read!(io, Lj)
        Lv = Vector{Float64}(undef, array_length)
        read!(io, Lv)
        rhs = Vector{Float64}(undef, N)
        read!(io, rhs)
        return N, Li, Lj, Lv, rhs
    end
end

"""
    read_system_of_equations_from_file(input_path::String)

Read system of equation files for use with pymumps solver. Matrix indices arrays Li and Lj are 1-based.

Returns:
- N: Int64 - Number of unknowns. Size of large matrix is NxN. For Stokes-continuity equation 
  N = (xnum-1)*(ynum-1)*3 and for heat equation N = xnum*ynum where xnum is number of nodes 
  of basic grid in x-direction and ynum is number of nodes of basic grid in y-direction.
- Li: Vector{Int64} - Large matrix 1-based row indices for non-zero matrix elements
- Lj: Vector{Int64} - Large matrix 1-based column indices for non-zero matrix elements
- Lv: Vector{Float64} - Non-zero matrix values of discretized system of equations
- R: Vector{Float64} - Right-hand-side vector of discretized equations
"""
function read_system_of_equations_from_file(input_path::String)
    file_name = "system_of_equations.jld"
    jld_file_path = joinpath(input_path, file_name)
    
    jldopen(jld_file_path, "r") do file
        matrix_info = file["matrix_info"]
        N = Int(matrix_info[1])
        
        Lv = file["matrix_nonzero_values"]
        Li = file["matrix_index_i"]
        Lj = file["matrix_index_j"]
        rhs = file["rhs_vector"]
        
        return N, Li, Lj, Lv, rhs
    end
end

"""
    read_solution_vector_file(input_path::String)

Read and then delete solution vector file.

Returns:
- S: Vector{Float64} - Solution vector. For Stokes-continuity equation 
  N = (xnum-1)*(ynum-1)*3 and for heat equation N = xnum*ynum where xnum is number of nodes 
  of basic grid in x-direction and ynum is number of nodes of basic grid in y-direction.
- mumps_solver_status: String - Status of the solver ("Success" or "Failure")
"""
function read_solution_vector_file(input_path::String)
    file_name = "solution_vector.jld"
    file_path = joinpath(input_path, file_name)
    
    if isfile(file_path)
        jldopen(file_path, "r") do file
            S = file["solution_vector"]
            rm(file_path)
            return S, "Success"
        end
    else
        println("Did not find solution file...")
        return zeros(1), "Failure"
    end
end

end # module 