""" Get material IDs for lithosphere layers.
"""
module GetBasicLithosphereIDs

import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import EarthBox.Markers.MarkerMaterials.GetMaterialIDs: get_material_id_dicts
import EarthBox.SurfaceProcesses.Sealevel.RelativeBaseLevel.ReferenceLithosphere: 
    get_lithosphere_layer_base_depths, LithosphereThicknesses

""" Layer base depths structure.

# Fields
- `y_base_upr_crust::Float64`: Upper crust base depth in meters
- `y_base_lwr_crust::Float64`: Lower crust base depth in meters
- `y_base_upr_mantle_lithosphere::Float64`: Upper mantle lithosphere base depth in meters
- `y_base_mid_mantle_lithosphere::Float64`: Mid mantle lithosphere base depth in meters
- `y_base_lwr_mantle_lithosphere::Float64`: Lower mantle lithosphere base depth in meters
"""
struct LayerBaseDepths
    y_base_upr_crust::Float64
    y_base_lwr_crust::Float64
    y_base_upr_mantle_lithosphere::Float64
    y_base_mid_mantle_lithosphere::Float64
    y_base_lwr_mantle_lithosphere::Float64
end

""" Lithosphere material IDs structure.

# Fields
- `upper_continental_crust::Int16`: Upper continental crust material ID
- `lower_continental_crust::Int16`: Lower continental crust material ID
- `mantle_lithosphere::Int16`: Mantle lithosphere material ID
- `asthenosphere::Int16`: Asthenosphere material ID
"""
struct LithosphereMaterialIDs
    upper_continental_crust::Int16
    lower_continental_crust::Int16
    mantle_lithosphere::Int16
    asthenosphere::Int16
end

""" Make lithosphere material IDs tuple.

# Arguments
- `materials::Materials`: Materials object

# Returns
- `LithosphereMaterialIDs`: Lithosphere material IDs structure
"""
function get_lithosphere_id_tuple(materials::Materials)::LithosphereMaterialIDs
    (id_upr_crust, id_lwr_crust, id_mantle_lithosphere, id_asthenosphere) = 
        get_layer_material_ids(materials)

    return LithosphereMaterialIDs(
        id_upr_crust,
        id_lwr_crust,
        id_mantle_lithosphere,
        id_asthenosphere
    )
end

""" Get material IDs for lithosphere layers.

# Arguments
- `materials::Materials`: Materials object

# Returns
- `Tuple{Int16, Int16, Int16, Int16}`: Tuple of material IDs for upper crust, lower crust, 
  mantle lithosphere, and asthenosphere
"""
function get_layer_material_ids(materials::Materials)::Tuple{Int16, Int16, Int16, Int16}
    materials_dict = materials.materials
    domain_matid_dict, _type_matid_dict = get_material_id_dicts(materials_dict)

    id_upr_crust = domain_matid_dict["UpperContinentalCrust"]
    id_lwr_crust = domain_matid_dict["LowerContinentalCrust"]
    
    if id_lwr_crust == -1 || id_upr_crust == -1
        id_upr_crust = domain_matid_dict["ContinentalCrust"]
        id_lwr_crust = domain_matid_dict["ContinentalCrust"]
    end
    
    if id_lwr_crust == -1 || id_upr_crust == -1
        throw(ArgumentError("Continental crust material ID not found."))
    end

    id_mantle_lithosphere = domain_matid_dict["MantleLithosphere"]
    if id_mantle_lithosphere == -1
        id_mantle_lithosphere = domain_matid_dict["UpperMantleLithosphere"]
    end
    if id_mantle_lithosphere == -1
        id_mantle_lithosphere = domain_matid_dict["MiddleMantleLithosphere"]
    end
    if id_mantle_lithosphere == -1
        id_mantle_lithosphere = domain_matid_dict["LowerMantleLithosphere"]
    end
    if id_mantle_lithosphere == -1
        throw(ArgumentError("Mantle lithosphere material ID not found."))
    end

    id_asthenosphere = domain_matid_dict["Asthenosphere"]
    if id_asthenosphere == -1
        id_asthenosphere = id_mantle_lithosphere
    end
    if id_asthenosphere == -1
        throw(ArgumentError("Asthenosphere material ID not found."))
    end

    return (id_upr_crust, id_lwr_crust, id_mantle_lithosphere, id_asthenosphere)
end

""" Get material id at a given y-depth.

# Arguments
- `layer_base_depths::LayerBaseDepths`: Base depths of lithosphere layers
- `lithosphere_material_ids::LithosphereMaterialIDs`: Material IDs for lithosphere layers
- `y_location_meters::Float64`: Depth in meters

# Returns
- `Int16`: Material ID at a given depth
"""
function get_material_id_y_depth(
    layer_base_depths::LayerBaseDepths,
    lithosphere_material_ids::LithosphereMaterialIDs,
    y_location_meters::Float64
)::Int16
    if y_location_meters <= layer_base_depths.y_base_upr_crust
        return lithosphere_material_ids.upper_continental_crust
    elseif y_location_meters <= layer_base_depths.y_base_lwr_crust
        return lithosphere_material_ids.lower_continental_crust
    elseif y_location_meters <= layer_base_depths.y_base_upr_mantle_lithosphere
        return lithosphere_material_ids.mantle_lithosphere
    elseif y_location_meters <= layer_base_depths.y_base_mid_mantle_lithosphere
        return lithosphere_material_ids.mantle_lithosphere
    elseif y_location_meters <= layer_base_depths.y_base_lwr_mantle_lithosphere
        return lithosphere_material_ids.mantle_lithosphere
    else
        return lithosphere_material_ids.asthenosphere
    end
end

""" Calculate material id array.

# Arguments
- `gridy::Vector{Float64}`: Y-coordinate grid
- `lithosphere_thicknesses::LithosphereThicknesses`: Lithosphere thicknesses
- `layer_material_ids::LithosphereMaterialIDs`: Layer material IDs

# Returns
- `Vector{Int16}`: Material ID array
"""
function get_material_id_grid(
    gridy::Vector{Float64},
    lithosphere_thicknesses::LithosphereThicknesses,
    layer_material_ids::LithosphereMaterialIDs
)::Vector{Int16}
    layer_base_depths = get_layer_base_depths(lithosphere_thicknesses)
    nnodes = length(gridy)
    matid_gridy = zeros(Int16, nnodes)
    
    for i in 1:nnodes
        y_location_meters = gridy[i]
        material_id = get_material_id_y_depth(
            layer_base_depths, layer_material_ids, y_location_meters
        )
        matid_gridy[i] = material_id
    end
    
    return matid_gridy
end

""" Get base depths of lithosphere layers.

# Arguments
- `lith_thicknesses::LithosphereThicknesses`: Lithosphere thicknesses

# Returns
- `LayerBaseDepths`: Layer base depths structure
"""
function get_layer_base_depths(lith_thicknesses::LithosphereThicknesses)::LayerBaseDepths
    (
        y_base_upr_crust,
        y_base_lwr_crust,
        y_base_upr_mantle_lithosphere,
        y_base_mid_mantle_lithosphere,
        y_base_lwr_mantle_lithosphere
    ) = get_lithosphere_layer_base_depths(lith_thicknesses)
    
    return LayerBaseDepths(
        y_base_upr_crust,
        y_base_lwr_crust,
        y_base_upr_mantle_lithosphere,
        y_base_mid_mantle_lithosphere,
        y_base_lwr_mantle_lithosphere
    )
end

end # module
