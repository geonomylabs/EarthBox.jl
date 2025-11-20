module LithosphericExtensionWeakFault

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: set_model_data!
import ..InitManager.InitStructs: LayerThickness
import ..InitManager.InitStructs: SevenLayerMaterialIDs
import ..InitManager.InitStructs: WeakFaultGeometry
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
    weak_fault = model.geometry.parameters.weak_fault
    weak_fault_geometry = WeakFaultGeometry(
        weak_fault.fault_dip_degrees.value,
        weak_fault.fault_thickness.value,
        weak_fault.x_initial_fault.value,
        weak_fault.fault_height.value
    )

    earth_layering = model.geometry.parameters.earth_layering
    layer_thickness = LayerThickness(
        model.geometry.parameters.sticky_air_geometry.thick_air.value,
        earth_layering.thick_upper_crust.value,
        earth_layering.thick_lower_crust.value,
        earth_layering.thick_upper_lith.value,
        earth_layering.thick_middle_lith.value,
        earth_layering.thick_lower_lith.value
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

    matid_weak_crustal_fault = matid_domains_dict["WeakCrustalFaultZone"]
    matid_weak_mantle_fault = matid_domains_dict["WeakMantleFaultZone"]

    marknum = model.markers.parameters.distribution.marknum.value
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array

    for imarker in 1:marknum
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        matid = define_material(
            x_marker, y_marker, layer_thickness,
            seven_layer_material_ids, weak_fault_geometry,
            matid_weak_crustal_fault, matid_weak_mantle_fault
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
    weak_fault_geometry::WeakFaultGeometry,
    matid_weak_crustal_fault::Int16,
    matid_weak_mantle_fault::Int16
)::Int16
    matid = SevenLayerEarthModel2D.define_material(
        y_marker,
        layer_thickness,
        seven_layer_material_ids
    )

    matid = define_weak_fault_region(
        matid,
        x_marker,
        y_marker,
        weak_fault_geometry,
        matid_weak_crustal_fault,
        matid_weak_mantle_fault,
        layer_thickness
    )

    return matid
end

function define_weak_fault_region(
    matid_ini::Int16,
    x_marker::Float64,
    y_marker::Float64,
    weak_fault_geometry::WeakFaultGeometry,
    matid_weak_crustal_fault::Int16,
    matid_weak_mantle_fault::Int16,
    layer_thickness::LayerThickness
)::Int16
    fault_dip_degrees = weak_fault_geometry.fault_dip_degrees
    fault_dip = fault_dip_degrees * π / 180.0
    dy_dx = tan(fault_dip)

    in_fault = in_weak_fault(
        x_marker,
        y_marker,
        dy_dx,
        weak_fault_geometry,
        layer_thickness
    )

    matid = matid_ini
    # weak zone (crust)
    if in_fault == 1
        matid = matid_weak_crustal_fault
    # weak zone (mantle)
    elseif in_fault == 2
        matid = matid_weak_mantle_fault
    end
    return matid
end

function in_weak_fault(
    x_marker::Float64,
    y_marker::Float64,
    dy_dx::Float64,
    weak_fault_geometry::WeakFaultGeometry,
    layer_thickness::LayerThickness
)::Int16
    vertical_fault_thickness = calculate_vertical_fault_thickness(weak_fault_geometry)
    y_moho = calculate_y_moho(layer_thickness)
    y_fault_limit = calculate_y_fault_limit(weak_fault_geometry, layer_thickness)

    x_initial_fault = weak_fault_geometry.x_initial_fault
    thick_air = layer_thickness.thick_air

    in_fault = 0
    if x_marker >= x_initial_fault && thick_air <= y_marker <= y_fault_limit
        dx_fault_base = x_marker - x_initial_fault
        y_fault_base = thick_air + dy_dx * dx_fault_base
        y_fault_top = y_fault_base - vertical_fault_thickness
        if y_fault_top <= y_marker <= y_fault_base
            in_fault = 1
            if y_marker > y_moho
                in_fault = 2
            end
        end
    end
    return in_fault
end

function calculate_y_fault_limit(
    weak_fault_geometry::WeakFaultGeometry,
    layer_thickness::LayerThickness
)::Float64
    thick_air = layer_thickness.thick_air
    fault_height = weak_fault_geometry.fault_height
    y_fault_limit = thick_air + fault_height
    return y_fault_limit
end

function calculate_vertical_fault_thickness(
    weak_fault_geometry::WeakFaultGeometry
)::Float64
    fault_thickness = weak_fault_geometry.fault_thickness
    fault_dip_degrees = weak_fault_geometry.fault_dip_degrees
    vertical_fault_thickness = fault_thickness * cos(fault_dip_degrees * π / 180.0)
    return vertical_fault_thickness
end

function calculate_y_moho(
    layer_thickness::LayerThickness
)::Float64
    thick_air = layer_thickness.thick_air
    thick_upper_crust = layer_thickness.thick_upper_crust
    thick_lower_crust = layer_thickness.thick_lower_crust
    y_moho = thick_air + thick_upper_crust + thick_lower_crust
    return y_moho
end

end # module