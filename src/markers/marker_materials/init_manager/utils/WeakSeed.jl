module WeakSeed

import ..InitManager.InitStructs: SeedCoordinates

"""
    define_material(
        matid_ini::Int,
        mxx::Float64,
        myy::Float64,
        seed_coordinates::SeedCoordinates,
        matid_weak_seed_mantle::Int
    )::Int

Define marker material ID for weak seed region.

# Arguments
- `matid_ini`: Initial material ID
- `mxx`: X-coordinate of marker
- `myy`: Y-coordinate of marker
- `seed_coordinates`: Seed region coordinates
- `matid_weak_seed_mantle`: Material ID for weak seed mantle

# Returns
- Material ID of marker
"""
function define_material(
    matid_ini::Int16,
    mxx::Float64,
    myy::Float64,
    seed_coordinates::SeedCoordinates,
    matid_weak_seed_mantle::Int16
)::Int16
    x_seed = seed_coordinates.x_seed
    y_seed = seed_coordinates.y_seed
    w_seed = seed_coordinates.w_seed
    
    xmin_seed, xmax_seed, ymin_seed, ymax_seed = get_seed_limits_simple_seed(
        x_seed, y_seed, w_seed)
    
    matid = matid_ini
    if ymin_seed >= myy <= ymax_seed
        if xmin_seed <= mxx <= xmax_seed
            matid = matid_weak_seed_mantle
        end
    end
    
    return matid
end

"""
    get_seed_limits_simple_seed(
        x_seed::Float64,
        y_seed::Float64,
        w_seed::Float64
    )::Tuple{Float64, Float64, Float64, Float64}

Calculate limits of weak seed.

# Arguments
- `x_seed`: X-coordinate of seed center
- `y_seed`: Y-coordinate of seed center
- `w_seed`: Width of seed

# Returns
- Tuple containing (xmin_seed, xmax_seed, ymin_seed, ymax_seed)
"""
function get_seed_limits_simple_seed(
    x_seed::Float64,
    y_seed::Float64,
    w_seed::Float64
)::Tuple{Float64, Float64, Float64, Float64}
    xmin_seed = x_seed - w_seed / 2.0
    xmax_seed = x_seed + w_seed / 2.0
    ymin_seed = y_seed - w_seed / 2.0
    ymax_seed = y_seed + w_seed / 2.0
    return xmin_seed, xmax_seed, ymin_seed, ymax_seed
end

end # module WeakSeed 