module TimeSteps

import ...XdmfParts: get_xdmf_start_for_collection_grid_and_domain
import ...XdmfParts: get_xdmf_end_for_collection_grid_and_domain

mutable struct VelocityXdmfTimeSteps
    output_dir::String
    xdmf_filepath::String
    xdmf_string::String

    function VelocityXdmfTimeSteps(output_dir::String)
        xdmf_filepath = joinpath(output_dir, "velocity.xdmf")
        xdmf_string = get_xdmf_start_for_collection_grid_and_domain("velocity")
        new(output_dir, xdmf_filepath, xdmf_string)
    end
end

function add_step_to_xdmf_string(steps::VelocityXdmfTimeSteps, string::String)::Nothing
    steps.xdmf_string *= string
    return nothing
end

function save_xdmf_string(steps::VelocityXdmfTimeSteps)::Nothing
    xdmf_string_updated = steps.xdmf_string * 
        get_xdmf_end_for_collection_grid_and_domain()
    open(steps.xdmf_filepath, "w") do io
        write(io, xdmf_string_updated)
    end
    return nothing
end

end # module
