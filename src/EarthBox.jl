"""
    EarthBox

Main module for the EarthBox project.

"""
module EarthBox

# Load startup configuration
include("__init__.jl")
include("earthbox_entry/imports.jl")
include("earthbox_entry/exports.jl")

const PDATA = get_eb_parameters()

include("earthbox_entry/eb_init_params.jl")
include("earthbox_entry/eb_state.jl")

function get_base_path()
    return @__DIR__
end

include("earthbox_entry/eb_run.jl")

end # module