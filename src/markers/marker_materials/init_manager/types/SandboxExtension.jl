module SandboxExtension

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!
import EarthBox.MobileWall: define_mobile_wall, define_plate_extension
import ..LayeredSand: define_layered_structure_in_sand

function initialize!(
    model::ModelData;
    parameters::Union{Dict{String, Any}, Nothing}=nothing,
    material_domain_ids::Union{Dict{String, Int}, Nothing}=nothing
)::Nothing
    set_model_data!(model, parameters, material_domain_ids)
    initialize_material_ids!(model)
    return nothing
end

function initialize_material_ids!(model::ModelData)::Nothing
    mobile_wall = model.geometry.parameters.mobile_wall
    x_left_mobile_wall = mobile_wall.x_left_mobile_wall.value
    x_right_mobile_wall = mobile_wall.x_right_mobile_wall.value
    y_top_mobile_wall = mobile_wall.y_top_mobile_wall.value
    y_bottom_mobile_wall = mobile_wall.y_bottom_mobile_wall.value
    plate_extension_thickness = mobile_wall.plate_extension_thickness.value
    plate_extension_width = mobile_wall.plate_extension_width.value

    sandbox = model.geometry.parameters.sandbox
    nsand_layers = sandbox.nsand_layers.value
    y_sand_air_interface = sandbox.y_sand_air_interface.value
    pdms_layer_width = sandbox.pdms_layer_width.value
    pdms_layer_thickness = sandbox.pdms_layer_thickness.value

    ysize = model.grids.parameters.geometry.ysize.value

    domains = model.materials.dicts.matid_domains
    matid_air = domains["Atmosphere"]
    matid_layered_sand_a = domains["SandA"]
    matid_layered_sand_b = domains["SandB"]
    matid_pdms_layer = domains["PDMSLayer"]
    matid_mobile_wall = domains["MobileWall"]
    matid_plate_extension = domains["PlateExtension"]

    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]

        matid = define_layered_structure_in_sand(
            matid_layered_sand_a, matid_layered_sand_b,
            y_marker, ysize, nsand_layers
        )

        matid = define_pdms_layer(
            matid, matid_pdms_layer, x_marker, y_marker, x_left_mobile_wall,
            y_bottom_mobile_wall, pdms_layer_width, pdms_layer_thickness
        )

        matid = define_mobile_wall(
            matid, matid_mobile_wall, x_marker, y_marker,
            x_left_mobile_wall, x_right_mobile_wall,
            y_top_mobile_wall, y_bottom_mobile_wall
        )

        matid = define_plate_extension(
            matid, x_marker, y_marker, matid_plate_extension,
            plate_extension_thickness, plate_extension_width,
            x_left_mobile_wall, y_bottom_mobile_wall
        )

        matid = define_air(
            matid, matid_air, x_marker, y_marker,
            y_sand_air_interface, x_left_mobile_wall, x_right_mobile_wall,
            y_top_mobile_wall, y_bottom_mobile_wall
        )
        model.markers.arrays.material.marker_matid.array[imarker] = matid
    end

    return nothing
end

function define_pdms_layer(
    matid::Int16,
    matid_pdms_layer::Int16,
    x_marker::Float64,
    y_marker::Float64,
    x_left_mobile_wall::Float64,
    y_bottom_pdms_layer::Float64,
    pdms_layer_width::Float64,
    pdms_layer_thickness::Float64
)::Int16
    xmid = x_left_mobile_wall/2.0
    xmin = xmid - pdms_layer_width/2.0
    xmax = xmid + pdms_layer_width/2.0
    ymin = y_bottom_pdms_layer - pdms_layer_thickness
    ymax = y_bottom_pdms_layer
    if xmin < x_marker < xmax && ymin < y_marker < ymax
        matid = matid_pdms_layer
    end
    return matid
end

function define_air(
    matid::Int16,
    matid_air::Int16,
    x_marker::Float64,
    y_marker::Float64,
    y_sand_air_interface::Float64,
    x_left_mobile_wall::Float64,
    x_right_mobile_wall::Float64,
    y_top_mobile_wall::Float64,
    y_bottom_mobile_wall::Float64
)::Int16
    if x_marker >= x_right_mobile_wall && y_marker <= y_bottom_mobile_wall
        matid = matid_air
    end
    if y_marker <= y_top_mobile_wall
        matid = matid_air
    end
    if x_marker < x_left_mobile_wall && y_marker <= y_sand_air_interface
        matid = matid_air
    end
    return matid
end

end # module SimpleExtension 