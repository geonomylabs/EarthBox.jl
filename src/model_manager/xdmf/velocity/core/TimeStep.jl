module TimeStep

import JLD2
import EarthBox.EarthBoxDtypes: AbstractXDMFTimeStep
import ...OutputDTypes: Mesh2djld, Vector2D
import ...XdmfParts: get_xdmf_start_of_timestep_grid
import ...XdmfParts: get_xdmf_end_of_timestep_grid
import ...XdmfParts: get_xdmf_topology_2drectmesh
import ...XdmfParts: get_xdmf_geometry_2d_vxvy
import ...XdmfParts: get_xdmf_vector_attribute_on_nodes
import ...XdmfTimeStepTools

struct VelocityXdmfTimeStep <: AbstractXDMFTimeStep
    model_time::Float64
    time_units::String
    noutput::Int
    mesh2djld::Mesh2djld
    vector::Vector2D
end

function get_xdmf_string_for_timestep(data::VelocityXdmfTimeStep)::String
    string = get_xdmf_start_of_timestep_grid(
        "velocity", data.model_time, data.time_units)
    string *= get_xdmf_topology_2drectmesh(
        data.mesh2djld.ynum, data.mesh2djld.xnum)
    string *= get_xdmf_geometry_2d_vxvy(
        data.mesh2djld.ynum, data.mesh2djld.xnum,
        data.mesh2djld.jld_mesh_filename, data.mesh2djld.jld_dataname_1dx,
        data.mesh2djld.jld_dataname_1dy
    )
    string *= get_xdmf_vector_attribute_on_nodes(
        data.mesh2djld.ynum,
        data.mesh2djld.xnum,
        data.mesh2djld.jld_mesh_filename,
        data.vector
    )
    string *= get_xdmf_end_of_timestep_grid()
    return string
end

function make_jld2_file(data::VelocityXdmfTimeStep, output_dir::String)::Nothing
    jld_mesh_filename = data.mesh2djld.jld_mesh_filename
    jld_mesh_filepath = joinpath(output_dir, jld_mesh_filename)
    
    JLD2.jldopen(jld_mesh_filepath, "w") do jldfile
        XdmfTimeStepTools.create_grid_entry_for_jld_file(data, jldfile)
        create_velocity_entry_for_jld_file(data, jldfile)
    end
    return nothing
end

function create_velocity_entry_for_jld_file(
    data::VelocityXdmfTimeStep,
    jldfile::JLD2.JLDFile
)::Nothing
    # Create velocity XYZ dataset
    group = JLD2.Group(jldfile, data.vector.jld_dataname_vxyz)
    group["array"] = data.vector.vectorxyz
    group["name"] = data.vector.name
    group["units"] = data.vector.units
    group["grid_type"] = data.vector.grid_type
    XdmfTimeStepTools.define_time_attr(data, group)

    # Create velocity X dataset
    group = JLD2.Group(jldfile, data.vector.jld_dataname_vx)
    group["array"] = data.vector.vectorx
    group["name"] = data.vector.name
    group["units"] = data.vector.units
    group["grid_type"] = data.vector.grid_type
    XdmfTimeStepTools.define_time_attr(data, group)

    # Create velocity Y dataset
    group = JLD2.Group(jldfile, data.vector.jld_dataname_vy)
    group["array"] = data.vector.vectory
    group["name"] = data.vector.name
    group["units"] = data.vector.units
    group["grid_type"] = data.vector.grid_type
    XdmfTimeStepTools.define_time_attr(data, group)

    # Create velocity magnitude dataset
    group = JLD2.Group(jldfile, data.vector.jld_dataname_vmag)
    group["array"] = data.vector.vmag
    group["name"] = data.vector.name
    group["units"] = data.vector.units
    group["grid_type"] = data.vector.grid_type
    XdmfTimeStepTools.define_time_attr(data, group)
    
    return nothing
end

end