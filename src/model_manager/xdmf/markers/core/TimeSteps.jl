module TimeSteps

import ...XdmfParts: get_xdmf_start_for_collection_grid_and_domain
import ...XdmfParts: get_xdmf_end_for_collection_grid_and_domain

mutable struct MarkersXdmfTimeSteps
    output_dir::String
    xdmf_filepath::String
    xdmf_string::String

    function MarkersXdmfTimeSteps(output_dir::String)
        xdmf_filepath = joinpath(output_dir, "markers.xdmf")
        xdmf_string = get_xdmf_start_for_collection_grid_and_domain("fields")
        new(output_dir, xdmf_filepath, xdmf_string)
    end
end

function add_step_to_xdmf_string(markers_xdmf::MarkersXdmfTimeSteps, string::String)
    markers_xdmf.xdmf_string *= string
end

function save(markers_xdmf::MarkersXdmfTimeSteps)
    xdmf_string_updated = markers_xdmf.xdmf_string * 
                         get_xdmf_end_for_collection_grid_and_domain()
    open(markers_xdmf.xdmf_filepath, "w") do io
        write(io, xdmf_string_updated)
    end
end

end