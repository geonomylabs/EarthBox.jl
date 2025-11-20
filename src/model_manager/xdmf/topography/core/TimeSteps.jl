module TimeSteps

import ...XdmfParts: get_xdmf_start_for_collection_grid_and_domain
import ...XdmfParts: get_xdmf_end_for_collection_grid_and_domain

mutable struct TopographyXdmfTimeSteps
    output_dir::String
    xdmf_filepath::String
    xdmf_string::String
    
    function TopographyXdmfTimeSteps(output_dir::String)
        xdmf_filepath = joinpath(output_dir, "topography.xdmf")
        xdmf_string = get_xdmf_start_for_collection_grid_and_domain("topo")
        new(output_dir, xdmf_filepath, xdmf_string)
    end
end

function add_step_to_xdmf_string(steps::TopographyXdmfTimeSteps, string::String)
    steps.xdmf_string *= string
end

function save(steps::TopographyXdmfTimeSteps)
    xdmf_string_updated = steps.xdmf_string * 
        get_xdmf_end_for_collection_grid_and_domain()
    open(steps.xdmf_filepath, "w") do io
        write(io, xdmf_string_updated)
    end
end

end # module
