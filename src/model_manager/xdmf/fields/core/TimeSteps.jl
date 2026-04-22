module TimeSteps

import ...XdmfParts: get_xdmf_start_for_collection_grid_and_domain
import ...XdmfParts: get_xdmf_end_for_collection_grid_and_domain


mutable struct FieldsXdmfTimeSteps
    output_dir::String
    xdmf_filepath::String
    xdmf_chunks::Vector{String}
end

function FieldsXdmfTimeSteps(output_dir::String)::FieldsXdmfTimeSteps
    xdmf_filepath = joinpath(output_dir, "fields.xdmf")
    chunks = [get_xdmf_start_for_collection_grid_and_domain("fields")]
    return FieldsXdmfTimeSteps(output_dir, xdmf_filepath, chunks)
end

function add_step_to_xdmf_string!(data::FieldsXdmfTimeSteps, string::String)::Nothing
    push!(data.xdmf_chunks, string)
    return nothing
end

function save(data::FieldsXdmfTimeSteps)::Nothing
    open(data.xdmf_filepath, "w") do io
        for chunk in data.xdmf_chunks
            write(io, chunk)
        end
        write(io, get_xdmf_end_for_collection_grid_and_domain())
    end
    return nothing
end

end
