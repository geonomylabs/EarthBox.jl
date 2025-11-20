module XdmfTimeStepTools

using JLD2
import JLD2: attributes, create_dataset
import EarthBox.EarthBoxDtypes: AbstractXDMFTimeStep
import ..OutputDTypes: Mesh2djld

function set_time(
    data::AbstractXDMFTimeStep,
    jldfile::JLD2.JLDFile
)
    jldfile["Myr"] = data.model_time
    jldfile["noutput"] = data.noutput
end

function set_y_sealevel(
    data::AbstractXDMFTimeStep,
    jldfile::JLD2.JLDFile
)
    jldfile["y_sealevel_meters"] = data.y_sealevel
end

function set_base_level_shift(
    data::AbstractXDMFTimeStep,
    jldfile::JLD2.JLDFile
)
    jldfile["base_level_shift_meters"] = data.base_level_shift
end

function define_time_attr(
    data::AbstractXDMFTimeStep,
    group::JLD2.Group
)
    group["noutput"] = data.noutput
    group["model_time"] = data.model_time
    group["time_units"] = data.time_units
end

""" Create basic and pressure grid entries for jld2 file.

Note that the sign of the y-axis is flipped for compatibility with Paraview.
"""
function create_grid_entry_for_jld_file(
    data::AbstractXDMFTimeStep,
    jldfile::JLD2.JLDFile
)::Nothing
    set_time(data, jldfile)

    jld_dataname_1dy = data.mesh2djld.jld_dataname_1dy
    group = JLD2.Group(jldfile, jld_dataname_1dy)
    group["array"] = -1.0 .* data.mesh2djld.gridy_km_array
    group["name"] = "gridy"
    group["units"] = "km"
    define_time_attr(data, group)

    jld_dataname_1dx = data.mesh2djld.jld_dataname_1dx
    group = JLD2.Group(jldfile, jld_dataname_1dx)
    group["array"] = data.mesh2djld.gridx_km_array
    group["name"] = "gridx"
    group["units"] = "km"
    define_time_attr(data, group)

    jld_dataname_pr1dx = data.mesh2djld.jld_dataname_pr1dx
    group = JLD2.Group(jldfile, jld_dataname_pr1dx)
    group["array"] = data.mesh2djld.gridx_pr_km_array
    group["name"] = "gridx_pr"
    group["units"] = "km"
    define_time_attr(data, group)

    jld_dataname_pr1dy = data.mesh2djld.jld_dataname_pr1dy
    group = JLD2.Group(jldfile, jld_dataname_pr1dy)
    group["array"] = -1.0 .* data.mesh2djld.gridy_pr_km_array
    group["name"] = "gridy_pr"
    group["units"] = "km"
    define_time_attr(data, group)

    return nothing
end

end # module
