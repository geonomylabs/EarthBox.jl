module SystemExporter

import EarthBox.PrintFuncs: print_info
using JLD2

function export_system_of_equations_binary(
    filepath::String, 
    N::Int, 
    Li::Vector{Int64}, 
    Lj::Vector{Int64}, 
    Lv::Vector{Float64}, 
    rhs::Vector{Float64}
)::Nothing
    open(filepath, "w") do io
        write(io, N)
        write(io, length(Li))
        write(io, Li)
        write(io, Lj)
        write(io, Lv)
        write(io, rhs)
    end
    return nothing
end

"""
    export_system_of_equations_for_mumps(
        output_dir::String,
        N::Int,
        Li::Vector{Int64},
        Lj::Vector{Int64},
        Lv::Vector{Float64},
        R::Vector{Float64}
    ) -> String

Save system of equations (soe) to jld2 files for use with pymumps.

# Arguments
- `output_dir`: Directory to save the output files
- `N`: Matrix dimension
- `Li`: Row indices of non-zero elements
- `Lj`: Column indices of non-zero elements
- `Lv`: Values of non-zero elements
- `R`: Right-hand side vector

# Returns
- Path to the system of equations directory
"""
function export_system_of_equations_for_mumps(
    output_dir::String,
    N::Int,
    Li::Vector{Int64},
    Lj::Vector{Int64},
    Lv::Vector{Float64},
    R::Vector{Float64},
    dir_name::Union{String, Nothing} = "system_of_equations"
)::String
    if dir_name === nothing
        soe_dir_path = output_dir
    else
        soe_dir_path = joinpath(output_dir, dir_name)
    end
    if !isdir(soe_dir_path)
        mkpath(soe_dir_path)
    end
    matrix_info = [N, N]
    jld_file_path = joinpath(soe_dir_path, "system_of_equations.jld")
    print_info("Writing system of equations to $jld_file_path", level=2)
    jldopen(jld_file_path, "w") do file
        file["matrix_info"] = matrix_info
        file["matrix_nonzero_values"] = Lv
        file["matrix_index_i"] = Li
        file["matrix_index_j"] = Lj
        file["rhs_vector"] = R
    end
    return soe_dir_path
end

"""
    export_solution_vector(soe_dir_path::String, S::Vector{Float64})

Export solution vector for external solver.

# Arguments
- `soe_dir_path`: Directory containing the system of equations
- `S`: Solution vector to export
"""
function export_solution_vector(soe_dir_path::String, S::Vector{Float64})
    jld_name = "solution_vector.jld"
    jld_file_path = joinpath(soe_dir_path, jld_name)
    
    jldopen(jld_file_path, "w") do file
        file["solution_vector"] = S
    end
end

end # module 