module Sticky

import EarthBox.ModelDataContainer: ModelData

function get_sticky_material_ids(model::ModelData)::Tuple{Int16, Int16}
    types = model.materials.dicts.matid_types
    matid_sticky_air = types["StickyAir"][1]
    matid_sticky_water = types["StickyWater"][1]
    return (matid_sticky_air, matid_sticky_water)
end

function check_recycle_sticky(
    mitype::Int16,
    matids_sticky::Tuple{Int16, Int16},
    nrecycle_sticky::Int
)::Bool
    insticky = check_insticky(mitype, matids_sticky)
    recycle_sticky = false
    if insticky && nrecycle_sticky > 0
        recycle_sticky = true
    end
    return recycle_sticky
end

function check_insticky(
    matid::Int16,
    matids_sticky::Tuple{Int16, Int16}
)::Bool
    insticky = false
    if matid in matids_sticky
        insticky = true
    end
    return insticky
end

end # module 