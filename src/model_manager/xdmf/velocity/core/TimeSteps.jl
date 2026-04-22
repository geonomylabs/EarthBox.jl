module TimeSteps

import ...XdmfParts: get_xdmf_start_for_collection_grid_and_domain
import ...XdmfParts: get_xdmf_end_for_collection_grid_and_domain

mutable struct VelocityXdmfTimeSteps
    output_dir::String
    xdmf_filepath::String
    xdmf_chunks::Vector{String}

    function VelocityXdmfTimeSteps(output_dir::String)
        xdmf_filepath = joinpath(output_dir, "velocity.xdmf")
        chunks = [get_xdmf_start_for_collection_grid_and_domain("velocity")]
        new(output_dir, xdmf_filepath, chunks)
    end
end

function add_step_to_xdmf_string(steps::VelocityXdmfTimeSteps, string::String)::Nothing
    push!(steps.xdmf_chunks, string)
    return nothing
end

function save_xdmf_string(steps::VelocityXdmfTimeSteps)::Nothing
    open(steps.xdmf_filepath, "w") do io
        for chunk in steps.xdmf_chunks
            write(io, chunk)
        end
        write(io, get_xdmf_end_for_collection_grid_and_domain())
    end
    return nothing
end

end # module
