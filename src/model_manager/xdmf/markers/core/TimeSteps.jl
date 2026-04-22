module TimeSteps

import ...XdmfParts: get_xdmf_start_for_collection_grid_and_domain
import ...XdmfParts: get_xdmf_end_for_collection_grid_and_domain

mutable struct MarkersXdmfTimeSteps
    output_dir::String
    xdmf_filepath::String
    xdmf_chunks::Vector{String}

    function MarkersXdmfTimeSteps(output_dir::String)
        xdmf_filepath = joinpath(output_dir, "markers.xdmf")
        chunks = [get_xdmf_start_for_collection_grid_and_domain("fields")]
        new(output_dir, xdmf_filepath, chunks)
    end
end

function add_step_to_xdmf_string(markers_xdmf::MarkersXdmfTimeSteps, string::String)
    push!(markers_xdmf.xdmf_chunks, string)
end

function save(markers_xdmf::MarkersXdmfTimeSteps)
    open(markers_xdmf.xdmf_filepath, "w") do io
        for chunk in markers_xdmf.xdmf_chunks
            write(io, chunk)
        end
        write(io, get_xdmf_end_for_collection_grid_and_domain())
    end
end

end