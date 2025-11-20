module SevenLayerEarthModel2D

import ..InitManager.InitStructs: LayerThickness, SevenLayerMaterialIDs
import ....Layers: calc_depths_for_layers

"""
    define_material(myy::Float64, layer_thickness::LayerThickness,
                   seven_layer_material_ids::SevenLayerMaterialIDs)::Int

Define material ID's for seven layer earth model.

# Arguments
- `myy`: Y-coordinate of marker
- `layer_thickness`: Layer thickness parameters
- `seven_layer_material_ids`: Material IDs for each layer

# Returns
- Material ID of marker
"""
function define_material(
    myy::Float64, 
    layer_thickness::LayerThickness,
    seven_layer_material_ids::SevenLayerMaterialIDs
)::Int16
    thick_layer1 = layer_thickness.thick_air
    thick_layer2 = layer_thickness.thick_upper_crust
    thick_layer3 = layer_thickness.thick_lower_crust
    thick_layer4 = layer_thickness.thick_upper_lith
    thick_layer5 = layer_thickness.thick_middle_lith
    thick_layer6 = layer_thickness.thick_lower_lith

    matid_layer1 = seven_layer_material_ids.matid_sticky_air
    matid_layer2 = seven_layer_material_ids.matid_upper_crust
    matid_layer3 = seven_layer_material_ids.matid_lower_crust
    matid_layer4 = seven_layer_material_ids.matid_upper_lith
    matid_layer5 = seven_layer_material_ids.matid_middle_lith
    matid_layer6 = seven_layer_material_ids.matid_lower_lith
    matid_layer7 = seven_layer_material_ids.matid_asthenosphere

    depth_layer1, depth_layer2, depth_layer3,
    depth_layer4, depth_layer5, depth_layer6 = calc_depths_for_layers(
        thick_layer1, thick_layer2,
        thick_layer3, thick_layer4,
        thick_layer5, thick_layer6
    )

    matid = matid_layer7
    if myy <= depth_layer1
        matid = matid_layer1
    elseif myy > depth_layer1 && myy <= depth_layer2
        matid = matid_layer2
    elseif myy > depth_layer2 && myy <= depth_layer3
        matid = matid_layer3
    elseif myy > depth_layer3 && myy <= depth_layer4
        matid = matid_layer4
    elseif myy > depth_layer4 && myy <= depth_layer5
        matid = matid_layer5
    elseif myy > depth_layer5 && myy <= depth_layer6
        matid = matid_layer6
    end
    
    return matid
end

end # module SevenLayerEarthModel2D 