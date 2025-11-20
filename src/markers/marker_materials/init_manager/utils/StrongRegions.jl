module StrongRegions

import ....Layers: calc_depths_for_layers
import ..InitManager.InitStructs: LayerThickness
import ..InitManager.InitStructs: StrongZoneLimits
import ..InitManager.InitStructs: StrongZoneMaterialIDs

"""
    define_material(
        x_marker::Float64,
        y_marker::Float64,
        matid_ini::Int,
        strong_zone_limits::StrongZoneLimits,
        strong_zone_material_ids::StrongZoneMaterialIDs,
        layer_thickness::LayerThickness
    )::Int

Define marker material ID's based on strong lithospheric zones.

# Arguments
- `x_marker`: X-coordinate of marker
- `y_marker`: Y-coordinate of marker
- `matid_ini`: Initial material ID
- `strong_zone_limits`: Strong zone boundary limits
- `strong_zone_material_ids`: Material IDs for strong zones
- `layer_thickness`: Layer thickness parameters

# Returns
- Material ID of marker
"""
function define_material(
    x_marker::Float64,
    y_marker::Float64,
    matid_ini::Int16,
    strong_zone_limits::StrongZoneLimits,
    strong_zone_material_ids::StrongZoneMaterialIDs,
    layer_thickness::LayerThickness
)::Int16
    thick_air = layer_thickness.thick_air
    thick_upper_crust = layer_thickness.thick_upper_crust
    thick_lower_crust = layer_thickness.thick_lower_crust
    thick_upper_lith = layer_thickness.thick_upper_lith
    thick_middle_lith = layer_thickness.thick_middle_lith
    thick_lower_lith = layer_thickness.thick_lower_lith

    (
        depth_air, depth_upper_crust, depth_lower_crust,
        depth_upper_lith, depth_middle_lith, depth_lower_lith
    ) = calc_depths_for_layers(
        thick_air, thick_upper_crust,
        thick_lower_crust, thick_upper_lith,
        thick_middle_lith, thick_lower_lith
    )

    x_left_strong = strong_zone_limits.x_left_strong
    x_right_strong = strong_zone_limits.x_right_strong

    matid = matid_ini
    if depth_air < y_marker <= depth_upper_crust
        if x_marker < x_left_strong || x_marker > x_right_strong
            matid = strong_zone_material_ids.matid_strong_upper_crust
        end
    elseif depth_upper_crust < y_marker <= depth_lower_crust
        if x_marker < x_left_strong || x_marker > x_right_strong
            matid = strong_zone_material_ids.matid_strong_lower_crust
        end
    elseif depth_lower_crust < y_marker <= depth_upper_lith
        if x_marker < x_left_strong || x_marker > x_right_strong
            matid = strong_zone_material_ids.matid_strong_lith
        end
    elseif depth_upper_lith < y_marker && y_marker <= depth_middle_lith
        if x_marker < x_left_strong || x_marker > x_right_strong
            matid = strong_zone_material_ids.matid_strong_lith
        end
    elseif depth_middle_lith < y_marker <= depth_lower_lith
        if x_marker < x_left_strong || x_marker > x_right_strong
            matid = strong_zone_material_ids.matid_strong_lith
        end
    end
    
    return matid
end

end # module StrongRegions 