module InitStructs

struct MaterialIniData
    xsize::Float64
    ysize::Float64
    thick_air::Float64
    thick_upper_crust::Float64
    thick_lower_crust::Float64
    thick_upper_lith::Float64
    thick_middle_lith::Float64
    thick_lower_lith::Float64
    x_seed::Float64
    y_seed::Float64
    w_seed::Float64
    depth_interface_h1::Float64
    xhole_start::Float64
    xhole_middle::Float64
    xhole_end::Float64
    xhole_depth::Float64
    fault_dip_degrees::Float64
    x_left_strong::Float64
    x_right_strong::Float64
    plume_radius::Float64
    plume_center_x::Float64
    plume_center_y::Float64
    plume_head_thick::Float64
end

struct LayerThickness
    thick_air::Float64
    thick_upper_crust::Float64
    thick_lower_crust::Float64
    thick_upper_lith::Float64
    thick_middle_lith::Float64
    thick_lower_lith::Float64
end

struct SeedCoordinates
    x_seed::Float64
    y_seed::Float64
    w_seed::Float64
end

struct TriangularHoleGeometry
    xhole_start::Float64
    xhole_middle::Float64
    xhole_end::Float64
    xhole_depth::Float64
end

struct StrongZoneLimits
    x_left_strong::Float64
    x_right_strong::Float64
end

struct PlumeGeometry
    plume_radius::Float64
    plume_center_x::Float64
    plume_center_y::Float64
    plume_head_thick::Float64
end

struct MaterialIDs
    matid_sticky_air::Int16
    matid_sticky_water::Int16
    matid_upper_crust::Int16
    matid_lower_crust::Int16
    matid_upper_lith::Int16
    matid_middle_lith::Int16
    matid_lower_lith::Int16
    matid_asthenosphere::Int16
    matid_weak_seed_crust::Int16
    matid_weak_seed_mantle::Int16
    matid_weak_crustal_fault::Int16
    matid_weak_mantle_fault::Int16
    matid_strong_upper_crust::Int16
    matid_strong_lower_crust::Int16
    matid_strong_lith::Int16
    matid_plume::Int16
end

struct SevenLayerMaterialIDs
    matid_sticky_air::Int16
    matid_upper_crust::Int16
    matid_lower_crust::Int16
    matid_upper_lith::Int16
    matid_middle_lith::Int16
    matid_lower_lith::Int16
    matid_asthenosphere::Int16
end

struct StrongZoneMaterialIDs
    matid_strong_upper_crust::Int
    matid_strong_lower_crust::Int
    matid_strong_lith::Int
end

struct WeakFaultGeometry
    fault_dip_degrees::Float64
    fault_thickness::Float64
    x_initial_fault::Float64
    fault_height::Float64
end

end # module InitializationStructs 