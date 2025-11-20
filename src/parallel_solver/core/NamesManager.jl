module NamesManager

"""
    FileAndDirNames()

Struct with file and directory names used in communication between the main process and the external 
MUMPS solver loop.

# Fields
- `termination_info_file_name`: Binary file with integer indicating whether external solver should terminate
    - 0 = don't terminate, -1 = terminate
- `ready_to_solve_file_name`: Text file signaling system is ready to be solved. The existence of this file
    indicates that the information has been exported necessary to solve the system of equations.
- `solution_flag_file_name`: Binary file with integer indicating solution status:
    - 0 = solution was not produced possibly due to an error in the solver loop, 
    - 1 = solution was successfully produced and exported
- `soe_dir_name`: Directory name where the system of equations and IO communicator files are stored
- `soe_file_name`: JLD file name containing the system of equations
"""
Base.@kwdef struct FileAndDirNames
    termination_info_file_name::String = "termination_info_0001.bin"
    ready_to_solve_file_name::String = "ready_to_solve_0001.txt"
    solution_flag_file_name::String = "solution_flag_0001.bin"
    soe_dir_name::String = "system_of_equations"
    soe_file_name::String = "system_of_equations.jld"
end

end # module