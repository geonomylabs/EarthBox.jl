module GetIDs

import EarthBox.ModelDataContainer: ModelData

function get_sticky_material_ids(model::ModelData)::Tuple{Int,Int}
    types = model.materials.dicts.matid_types
    matid_sticky_air = types["StickyAir"][1]
    matid_sticky_water = types["StickyWater"][1]
    return (matid_sticky_air, matid_sticky_water)
end 

end