module TimeSteps

import ...XdmfParts: get_xdmf_start_for_collection_grid_and_domain
import ...XdmfParts: get_xdmf_end_for_collection_grid_and_domain

mutable struct TopographyXdmfTimeSteps
    output_dir::String
    xdmf_filepath::String
    xdmf_chunks::Vector{String}

    function TopographyXdmfTimeSteps(output_dir::String)
        xdmf_filepath = joinpath(output_dir, "topography.xdmf")
        chunks = [get_xdmf_start_for_collection_grid_and_domain("topo")]
        new(output_dir, xdmf_filepath, chunks)
    end
end

function add_step_to_xdmf_string(steps::TopographyXdmfTimeSteps, string::String)
    push!(steps.xdmf_chunks, string)
end

function save(steps::TopographyXdmfTimeSteps)
    open(steps.xdmf_filepath, "w") do io
        for chunk in steps.xdmf_chunks
            write(io, chunk)
        end
        write(io, get_xdmf_end_for_collection_grid_and_domain())
    end
end

end # module
