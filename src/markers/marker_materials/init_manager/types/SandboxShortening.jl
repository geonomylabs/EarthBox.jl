module SandboxShortening

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!
import EarthBox.MobileWall: define_mobile_wall
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

    sandbox = model.geometry.parameters.sandbox
    nsand_layers = sandbox.nsand_layers.value
    y_sand_air_interface = sandbox.y_sand_air_interface.value
    y_top_microbeads = sandbox.y_top_microbeads.value
    y_bottom_microbeads = sandbox.y_bottom_microbeads.value
    x_left_ramp = sandbox.x_left_ramp.value
    y_left_ramp = y_sand_air_interface
    x_right_ramp = sandbox.x_right_ramp.value
    ramp_dip_deg = sandbox.ramp_dip_deg.value

    ysize = model.grids.parameters.geometry.ysize.value

    domains = model.materials.dicts.matid_domains
    matid_air = domains["Atmosphere"]
    matid_layered_sand_a = domains["SandA"]
    matid_layered_sand_b = domains["SandB"]
    matid_microbeads = domains["Microbeads"]
    matid_mobile_wall = domains["MobileWall"]

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

        matid = define_microbead_layer(
            matid, matid_microbeads,
            y_marker, y_top_microbeads, y_bottom_microbeads
        )

        matid = define_mobile_wall(
            matid, matid_mobile_wall, x_marker, y_marker,
            x_left_mobile_wall, x_right_mobile_wall,
            y_top_mobile_wall, y_bottom_mobile_wall
        )

        matid = define_air(
            matid, matid_air, x_marker, y_marker,
            x_left_ramp, x_right_mobile_wall, x_right_ramp,
            y_bottom_mobile_wall, y_left_ramp, y_top_mobile_wall,
            ramp_dip_deg
        )
        model.markers.arrays.material.marker_matid.array[imarker] = matid
    end

    return nothing
end

function define_microbead_layer(
    matid::Int16,
    matid_microbeads::Int16,
    y_marker::Float64,
    y_top_microbeads::Float64,
    y_bottom_microbeads::Float64
)::Int16
    if y_top_microbeads < y_marker < y_bottom_microbeads
        matid = matid_microbeads
    end
    return matid
end

function define_air(
    matid::Int16,
    matid_air::Int16,
    x_marker::Float64,
    y_marker::Float64,
    x_left_ramp::Float64,
    x_right_mobile_wall::Float64,
    x_right_ramp::Float64,
    y_bottom_mobile_wall::Float64,
    y_left_ramp::Float64,
    y_top_mobile_wall::Float64,
    ramp_dip_deg::Float64
)::Int16
    y_air_sand_interface = calculate_interface_depth(
        x_left_ramp, y_left_ramp, ramp_dip_deg, x_marker)

    if x_marker <= x_left_ramp && y_marker <= y_left_ramp
        matid = matid_air
    end

    if x_marker >= x_right_mobile_wall && y_marker <= y_bottom_mobile_wall
        matid = matid_air
    end

    if y_marker <= y_top_mobile_wall
        matid = matid_air
    end

    if x_left_ramp < x_marker <= x_right_ramp && y_marker <= y_air_sand_interface
        matid = matid_air
    end

    return matid
end

function calculate_interface_depth(
    x_left_ramp::Float64,
    y_left_ramp::Float64,
    ramp_dip_deg::Float64,
    x_marker::Float64
)::Float64
    return y_left_ramp - (x_marker - x_left_ramp)*tan(ramp_dip_deg/180.0*π)
end

end # module SandBoxShortening 