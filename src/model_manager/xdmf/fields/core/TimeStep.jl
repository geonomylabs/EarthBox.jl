module TimeStep

import JLD2
import EarthBox.EarthBoxDtypes: AbstractXDMFTimeStep
import ...OutputDTypes: Mesh2djld
import ...OutputDTypes: ScalarField
import ...XdmfParts: get_xdmf_start_of_timestep_grid
import ...XdmfParts: get_xdmf_end_of_timestep_grid
import ...XdmfParts: get_xdmf_topology_2drectmesh
import ...XdmfParts: get_xdmf_geometry_2d_vxvy
import ...XdmfParts: get_xdmf_scalar_attribute_on_nodes
import ...XdmfParts: get_xdmf_scalar_attribute_on_cells
import ...XdmfTimeStepTools

struct FieldsXdmfTimeStep <: AbstractXDMFTimeStep
    model_time::Float64
    time_units::String
    noutput::Int64
    mesh2djld::Mesh2djld
    scalars_on_nodes::Vector{ScalarField}
    scalars_on_cells::Vector{ScalarField}
end

function get_xdmf_string_for_timestep(data::FieldsXdmfTimeStep)::String
    string = get_xdmf_start_of_timestep_grid(
        "fields", data.model_time, data.time_units)
    string *= get_xdmf_topology_2drectmesh(
        data.mesh2djld.ynum, data.mesh2djld.xnum)
    string *= get_xdmf_geometry_2d_vxvy(
        data.mesh2djld.ynum, data.mesh2djld.xnum,
        data.mesh2djld.jld_mesh_filename, data.mesh2djld.jld_dataname_1dx,
        data.mesh2djld.jld_dataname_1dy
    )
    
    for scalar_field in data.scalars_on_nodes
        string *= get_xdmf_scalar_attribute_on_nodes(
            data.mesh2djld.ynum, data.mesh2djld.xnum,
            data.mesh2djld.jld_mesh_filename, scalar_field
        )
    end
    
    for scalar_field in data.scalars_on_cells
        string *= get_xdmf_scalar_attribute_on_cells(
            data.mesh2djld.ynum, data.mesh2djld.xnum,
            data.mesh2djld.jld_mesh_filename, scalar_field
        )
    end
    
    string *= get_xdmf_end_of_timestep_grid()
    return string
end

function make_jld2_file(
    data::FieldsXdmfTimeStep,
    output_dir::String
)::Nothing
    jld_mesh_filename = data.mesh2djld.jld_mesh_filename
    jld_mesh_filepath = joinpath(output_dir, jld_mesh_filename)
    
    JLD2.jldopen(jld_mesh_filepath, "w") do jldfile
        XdmfTimeStepTools.create_grid_entry_for_jld_file(data, jldfile)
        for scalar_field in data.scalars_on_nodes
            create_scalar_field_entry_for_jld_file(data, jldfile, scalar_field)
        end
        for scalar_field in data.scalars_on_cells
            create_scalar_field_entry_for_jld_file(data, jldfile, scalar_field)
        end
    end
end

function create_scalar_field_entry_for_jld_file(
    data::FieldsXdmfTimeStep,
    jldfile::JLD2.JLDFile,
    scalar_field::ScalarField
)::Nothing
    group = JLD2.Group(jldfile, scalar_field.jld_dataname)
    group["array"] = scalar_field.scalar_array
    group["name"] = scalar_field.name
    group["units"] = scalar_field.units
    group["grid_type"] = scalar_field.grid_type
    XdmfTimeStepTools.define_time_attr(data, group)
    return nothing
end

end