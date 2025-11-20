module FlexureTriangularHole

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!
import ..InitManager.InitStructs: LayerThickness
import ..InitManager.InitStructs: SevenLayerMaterialIDs
import ..InitManager.InitStructs: TriangularHoleGeometry
import ..SevenLayerEarthModel2D

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
    earth_layering = model.geometry.parameters.earth_layering
    crustal_hole = model.geometry.parameters.crustal_hole

    thick_air = model.geometry.parameters.sticky_air_geometry.thick_air.value

    layer_thickness = LayerThickness(
        model.geometry.parameters.sticky_air_geometry.thick_air.value,
        earth_layering.thick_upper_crust.value,
        earth_layering.thick_lower_crust.value,
        earth_layering.thick_upper_lith.value,
        earth_layering.thick_middle_lith.value,
        earth_layering.thick_lower_lith.value
    )

    triangular_hole_geometry = TriangularHoleGeometry(
        crustal_hole.xhole_start.value,
        crustal_hole.xhole_middle.value,
        crustal_hole.xhole_end.value,
        crustal_hole.xhole_depth.value
    )

    matid_domains_dict = model.materials.dicts.matid_domains

    seven_layer_material_ids = SevenLayerMaterialIDs(
        matid_domains_dict["Atmosphere"],
        matid_domains_dict["UpperContinentalCrust"],
        matid_domains_dict["LowerContinentalCrust"],
        matid_domains_dict["UpperMantleLithosphere"],
        matid_domains_dict["MiddleMantleLithosphere"],
        matid_domains_dict["LowerMantleLithosphere"],
        matid_domains_dict["Asthenosphere"]
    )

    matid_sticky_water = matid_domains_dict["Ocean"]

    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        matid = define_material(
            x_marker, y_marker, layer_thickness,
            seven_layer_material_ids, matid_sticky_water,
            triangular_hole_geometry, thick_air
        )

        model.markers.arrays.material.marker_matid.array[imarker] = matid
    end

    return nothing
end

function define_material(
    x_marker::Float64,
    y_marker::Float64,
    layer_thickness::LayerThickness,
    seven_layer_material_ids::SevenLayerMaterialIDs,
    matid_sticky_water::Int16,
    triangular_hole_geometry::TriangularHoleGeometry,
    thick_water::Float64
)::Int16
    matid = SevenLayerEarthModel2D.define_material(
        y_marker,
        layer_thickness,
        seven_layer_material_ids
    )

    matid = define_hole_in_crust(
        matid,
        x_marker,
        y_marker,
        matid_sticky_water,
        thick_water,
        triangular_hole_geometry
    )

    return matid
end

function define_hole_in_crust(
    matid_ini::Int16,
    marker_x::Float64,
    marker_y::Float64,
    matid_sticky_water::Int16,
    thick_water::Float64,
    triangular_hole_geometry::TriangularHoleGeometry
)::Int16
    xhole_start = triangular_hole_geometry.xhole_start
    xhole_middle = triangular_hole_geometry.xhole_middle
    xhole_end = triangular_hole_geometry.xhole_end
    xhole_depth = triangular_hole_geometry.xhole_depth

    matid = matid_ini
    on_left_side = on_left_side_of_hole(
        marker_x, marker_y, xhole_start, xhole_middle, thick_water)
    if on_left_side
        depth = calc_depth_hole_left(
            marker_x, thick_water, xhole_depth, xhole_start, xhole_middle)
        if marker_y < depth
            matid = matid_sticky_water
        end
    end

    on_right_side = on_right_side_of_hole(
        marker_x, marker_y, xhole_middle, xhole_end, thick_water)
    if on_right_side
        depth = calc_depth_hole_right(
            marker_x, thick_water, xhole_depth, xhole_middle, xhole_end)
        if marker_y < depth
            matid = matid_sticky_water
        end
    end

    return matid
end

function on_left_side_of_hole(
    marker_x::Float64,
    marker_y::Float64,
    xhole_start::Float64,
    xhole_middle::Float64,
    thick_water::Float64
)::Bool
    check = false
    if xhole_start < marker_x <= xhole_middle && marker_y >= thick_water
        check = true
    end
    return check
end

function on_right_side_of_hole(
    marker_x::Float64,
    marker_y::Float64,
    xhole_middle::Float64,
    xhole_end::Float64,
    thick_water::Float64
)::Bool
    check = false
    if xhole_middle < marker_x <= xhole_end && marker_y >= thick_water
        check = true
    end
    return check
end

function calc_depth_hole_left(
    marker_x::Float64,
    thick_water::Float64,
    xhole_depth::Float64,
    xhole_start::Float64,
    xhole_middle::Float64
)::Float64
    depth = (
        thick_water
      + xhole_depth/(xhole_middle - xhole_start)*(marker_x - xhole_start)
    )
    return depth
end

function calc_depth_hole_right(
    marker_x::Float64,
    thick_water::Float64,
    xhole_depth::Float64,
    xhole_middle::Float64,
    xhole_end::Float64
)::Float64
    depth = (
        thick_water
      + (
            xhole_depth
          - xhole_depth/(xhole_end - xhole_middle)
              *(marker_x - xhole_middle)
          )
    )
    return depth
end

end # module FlexureTriangularHole 