module TimeSteps

import ...XdmfParts: get_xdmf_start_for_collection_grid_and_domain
import ...XdmfParts: get_xdmf_end_for_collection_grid_and_domain


mutable struct FieldsXdmfTimeSteps
    output_dir::String
    xdmf_filepath::String
    xdmf_string::String
end

function FieldsXdmfTimeSteps(output_dir::String)::FieldsXdmfTimeSteps
    xdmf_filepath = joinpath(output_dir, "fields.xdmf")
    xdmf_string = get_xdmf_start_for_collection_grid_and_domain("fields")
    return FieldsXdmfTimeSteps(output_dir, xdmf_filepath, xdmf_string)
end

function add_step_to_xdmf_string!(data::FieldsXdmfTimeSteps, string::String)::Nothing
    data.xdmf_string *= string
    return nothing
end

function save(data::FieldsXdmfTimeSteps)::Nothing
    xdmf_string_updated = data.xdmf_string * 
        get_xdmf_end_for_collection_grid_and_domain()
    open(data.xdmf_filepath, "w") do io
        write(io, xdmf_string_updated)
    end
    return nothing
end

end
