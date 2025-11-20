module WeakNotch

import ..InitManager.InitStructs: SeedCoordinates

""" Define marker ID for weak seed region.

# Arguments
- `matid_ini`: Initial material ID
- `mxx`: X-coordinate of marker
- `myy`: Y-coordinate of marker
- `depth_lower_crust`: Depth of lower crust
- `matid_weak_seed_mantle`: Material ID for weak seed mantle
- `seed_coordinates`: Seed region coordinates

# Returns
- Material ID of marker
"""
function define_material(
    matid_ini::Int16,
    mxx::Float64, myy::Float64,
    depth_lower_crust::Float64,
    matid_weak_seed_mantle::Int16,
    seed_coordinates::SeedCoordinates
)::Int16
    xmin_seed, xmax_seed, ymin_seed, ymax_seed = get_seed_limits_weak_notch(
        depth_lower_crust, seed_coordinates.x_seed, seed_coordinates.w_seed)
    
    if ymin_seed <= myy <= ymax_seed
        if xmin_seed <= mxx <= xmax_seed
            return matid_weak_seed_mantle
        end
    end
    return matid_ini
end

""" Get limits of weak seed.

# Arguments
- `depth_lower_crust`: Depth of lower crust
- `x_seed`: X-coordinate of seed center
- `w_seed`: Width of seed

# Returns
- Tuple containing (xmin_seed, xmax_seed, ymin_seed, ymax_seed)
"""
function get_seed_limits_weak_notch(
    depth_lower_crust::Float64,
    x_seed::Float64,
    w_seed::Float64
)::Tuple{Float64, Float64, Float64, Float64}
    xmin_seed = x_seed - w_seed / 2.0
    xmax_seed = x_seed + w_seed / 2.0
    ymin_seed = depth_lower_crust - w_seed / 2.0
    ymax_seed = depth_lower_crust

    (
        xmin_seed, xmax_seed, ymin_seed, ymax_seed
    ) = calculate_seed_limits_for_weak_notch(
        x_seed, depth_lower_crust, w_seed)

    return xmin_seed, xmax_seed, ymin_seed, ymax_seed
end

""" Calculate limits of weak seed.

# Arguments
- `x_seed`: X-coordinate of seed center
- `y_seed_base`: Base Y-coordinate of seed
- `w_seed`: Width of seed

# Returns
- Tuple containing (xmin_seed, xmax_seed, ymin_seed, ymax_seed)
"""
function calculate_seed_limits_for_weak_notch(
    x_seed::Float64,
    y_seed_base::Float64,
    w_seed::Float64
)::Tuple{Float64, Float64, Float64, Float64}
    xmin_seed = x_seed - w_seed / 2.0
    xmax_seed = x_seed + w_seed / 2.0
    ymin_seed = y_seed_base - w_seed / 2.0
    ymax_seed = y_seed_base
    return xmin_seed, xmax_seed, ymin_seed, ymax_seed
end

end # module