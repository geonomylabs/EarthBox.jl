module XdmfParts

import ..OutputDTypes: ScalarField, Vector2D

function get_xdmf_geometry_2d_vxvy(
    ynum::Int,
    xnum::Int,
    jld_mesh_filename::String,
    jld_dataname_1dx::String,
    jld_dataname_1dy::String
)::String
    tab_str1 = "\t"
    tab_str2 = "\t\t"
    string = "$(tab_str1)<Geometry GeometryType=\"VXVY\">\n"
    string *= "$(tab_str2)<DataItem Dimensions=\"$(xnum)\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"
    string *= "$(tab_str2)$(jld_mesh_filename):/$(jld_dataname_1dx)\n"
    string *= "$(tab_str2)</DataItem>\n"
    string *= "$(tab_str2)<DataItem Dimensions=\"$(ynum)\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"
    string *= "$(tab_str2)$(jld_mesh_filename):/$(jld_dataname_1dy)\n"
    string *= "$(tab_str2)</DataItem>\n"
    string *= "$(tab_str1)</Geometry>\n"
    return string
end

function get_xdmf_geometry_2d_xy_markers(
    nmarkers::Int,
    jld_markerfile::String,
    jld_dataname_xy::String
)::String
    tab_str1 = "\t"
    tab_str2 = "\t\t"
    string = "$(tab_str1)<Geometry GeometryType=\"XY\">\n"
    string *= "$(tab_str2)<DataItem Dimensions=\"$(nmarkers) 2\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n"
    string *= "$(tab_str2)$(jld_markerfile):/$(jld_dataname_xy)\n"
    string *= "$(tab_str2)</DataItem>\n"
    string *= "$(tab_str1)</Geometry>\n"
    return string
end

function get_xdmf_topology_2drectmesh(ynum::Int, xnum::Int)::String
    string = "\t<Topology TopologyType=\"2DRectMesh\" NumberOfElements=\"$(ynum) $(xnum)\"/>\n"
    return string
end

function get_xdmf_topology_polyvertex(
    nmarkers::Int,
    jld_markerfile::String,
    jld_dataname_ids::String
)::String
    tab_str1 = "\t"
    tab_str2 = "\t\t"
    string = "$(tab_str1)<Topology TopologyType=\"POLYVERTEX\" Dimensions=\"$(nmarkers)\" NumberOfElements=\"1\">\n"
    string *= "$(tab_str2)<DataItem Dimensions=\"$(nmarkers)\" NumberType=\"Int\" Format=\"HDF\">\n"
    string *= "$(tab_str2)$(jld_markerfile):/$(jld_dataname_ids)\n"
    string *= "$(tab_str2)</DataItem>\n"
    string *= "$(tab_str1)</Topology>\n"
    return string
end

function get_xdmf_topology_polyline(
    nmarkers::Int,
    jld_markerfile::String,
    jld_dataname_ids::String
)::String
    tab_str1 = "\t"
    string = "$(tab_str1)<Topology TopologyType=\"POLYLINE\" NodesPerElement=\"$(nmarkers)\">\n"
    string *= "$(tab_str1)</Topology>\n"
    return string
end

function get_xdmf_scalar_attribute_on_cells(
    ynum::Int,
    xnum::Int,
    jld_mesh_filename::String,
    scalar_field::ScalarField
)::String
    scalar_name = scalar_field.name
    jld_dataname_scalar = scalar_field.jld_dataname
    number_type = scalar_field.number_type
    tab_str1 = "\t"
    tab_str2 = "\t\t"
    string = "$(tab_str1)<Attribute Name=\"$(scalar_name)\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
    string *= "$(tab_str2)<DataItem Dimensions=\"$(ynum-1) $(xnum-1)\" NumberType=\"$(number_type)\" Precision=\"8\" Format=\"HDF\">\n"
    string *= "$(tab_str2)$(jld_mesh_filename):/$(jld_dataname_scalar)\n"
    string *= "$(tab_str2)</DataItem>\n"
    string *= "$(tab_str1)</Attribute>\n"
    return string
end

function get_xdmf_scalar_attribute_on_nodes(
    ynum::Int,
    xnum::Int,
    jld_mesh_filename::String,
    scalar_field::ScalarField
)::String
    scalar_name = scalar_field.name
    jld_dataname_scalar = scalar_field.jld_dataname
    number_type = scalar_field.number_type
    tab_str1 = "\t"
    tab_str2 = "\t\t"
    string = "$(tab_str1)<Attribute Name=\"$(scalar_name)\" AttributeType=\"Scalar\" Center=\"Node\">\n"
    string *= "$(tab_str2)<DataItem Dimensions=\"$(ynum) $(xnum)\" NumberType=\"$(number_type)\" Precision=\"8\" Format=\"HDF\">\n"
    string *= "$(tab_str2)$(jld_mesh_filename):/$(jld_dataname_scalar)\n"
    string *= "$(tab_str2)</DataItem>\n"
    string *= "$(tab_str1)</Attribute>\n"
    return string
end

function get_xdmf_scalar_attribute_on_nodes_for_markers(
    nmarkers::Int,
    jld_markerfile::String,
    scalar_field::ScalarField
)::String
    field_name = scalar_field.name
    jld_dataname_scalar = scalar_field.jld_dataname
    number_type = scalar_field.number_type
    tab_str1 = "\t"
    tab_str2 = "\t\t"
    string = "$(tab_str1)<Attribute Name=\"$(field_name)\" AttributeType=\"Scalar\" Center=\"Node\">\n"
    string *= "$(tab_str2)<DataItem Dimensions=\"$(nmarkers)\" NumberType=\"$(number_type)\" Precision=\"8\" Format=\"HDF\">\n"
    string *= "$(tab_str2)$(jld_markerfile):/$(jld_dataname_scalar)\n"
    string *= "$(tab_str2)</DataItem>\n"
    string *= "$(tab_str1)</Attribute>\n"
    return string
end

function get_xdmf_vector_attribute_on_nodes(
    ynum::Int,
    xnum::Int,
    jld_mesh_filename::String,
    vector::Vector2D
)::String
    vector_name = vector.name
    jld_dataname_vxyz = vector.jld_dataname_vxyz
    number_type = vector.number_type
    tab_str1 = "\t"
    tab_str2 = "\t\t"
    string = "$(tab_str1)<Attribute Name=\"$(vector_name)\" AttributeType=\"Vector\" Center=\"Node\">\n"
    string *= "$(tab_str2)<DataItem Dimensions=\"1 $(ynum) $(xnum) 3\" NumberType=\"$(number_type)\" Precision=\"8\" Format=\"HDF\">\n"
    string *= "$(tab_str2)$(jld_mesh_filename):/$(jld_dataname_vxyz)\n"
    string *= "$(tab_str2)</DataItem>\n"
    string *= "$(tab_str1)</Attribute>\n"
    return string
end

function get_xdmf_start_for_collection_grid_and_domain(grid_name::String)::String
    string = "<?xml version=\"1.0\" ?>\n"
    string *= "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n"
    string *= "<Domain>\n"
    string *= "<Grid Name=\"$(grid_name)\" GridType=\"Collection\" CollectionType=\"Temporal\">\n"
    return string
end

function get_xdmf_end_for_collection_grid_and_domain()::String
    string = "</Grid>\n"
    string *= "</Domain>\n"
    string *= "</Xdmf>\n"
    return string
end

function get_xdmf_start_of_timestep_grid(
    grid_name::String,
    model_time::Float64,
    units::String
)::String
    string = "\t<Grid Name=\"$(grid_name)\" GridType=\"Uniform\">\n"
    string *= "\t<Time Type=\"Single\" Units=\"$(units)\" Value=\"$(model_time)\"/>\n"
    return string
end

function get_xdmf_end_of_timestep_grid()::String
    string = "\t</Grid>\n"
    return string
end

function get_xdmf_start_for_grid_and_domain(grid_name::String)::String
    string = "<?xml version=\"1.0\" ?>\n"
    string *= "<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n"
    string *= "<Domain>\n"
    string *= "<Grid Name=\"$(grid_name)\" GridType=\"Uniform\">\n"
    return string
end

function get_xdmf_time_part(model_time::Float64)::String
    string = "\n\t<Time Type=\"Single\" Value=\"$(model_time)\" />\n\n"
    return string
end

function get_xdmf_end_for_grid_and_domain()::String
    string = "</Grid>\n"
    string *= "</Domain>\n"
    string *= "</Xdmf>\n"
    return string
end

end # module 