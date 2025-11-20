module DensityProps

include("DomainMaterialIDs.jl")

import .DomainMaterialIDs
import EarthBox.ModelDataContainer: ModelData
import EarthBox.Markers.MarkerMaterials: GetProperty
import EarthBox.Markers.MarkerTemperature.InitManager.TempIniStructs: LayerThickness

struct DensityProperties
    standard_density::Float64  # kg/m³
    expansivity::Float64      # 1/K
    compressibility::Float64  # 1/Pa
end

struct DensityModel
    upr_cont_crust::DensityProperties
    lwr_cont_crust::DensityProperties
    upr_mantle_lithosphere::DensityProperties
    mid_mantle_lithosphere::DensityProperties
    lwr_mantle_lithosphere::DensityProperties
    asthenosphere::DensityProperties
end

function get_density_model(model::ModelData)::DensityModel
    density_props_upr_cont_crust, density_props_lwr_cont_crust = 
        get_density_properties_continental_crust(model)

    density_props_upper_mantle_lithosphere, density_props_middle_mantle_lithosphere,
    density_props_lower_mantle_lithosphere = get_density_properties_mantle_lithosphere(model)

    density_props_asthenosphere = get_density_properties_asthenosphere(model)

    return DensityModel(
        density_props_upr_cont_crust,
        density_props_lwr_cont_crust,
        density_props_upper_mantle_lithosphere,
        density_props_middle_mantle_lithosphere,
        density_props_lower_mantle_lithosphere,
        density_props_asthenosphere
    )
end

function get_density_properties_asthenosphere(model::ModelData)::DensityProperties
    id_asthenosphere = DomainMaterialIDs.get_asthenosphere_material_id(model)
    if id_asthenosphere < 0
        return DensityProperties(3300.0, 3.0e-5, 1.0e-11)
    else
        return DensityProperties(
            GetProperty.get_standard_density(model, id_asthenosphere),
            GetProperty.get_thermal_expansivity(model, id_asthenosphere),
            GetProperty.get_compressibility(model, id_asthenosphere)
        )
    end
end

function get_density_properties_mantle_lithosphere(
    model::ModelData
)::Tuple{DensityProperties,DensityProperties,DensityProperties}
    id_upr_mantle_lithosphere, id_mid_mantle_lithosphere, id_lwr_mantle_lithosphere = 
        DomainMaterialIDs.get_lithosphere_material_ids(model)

    if id_upr_mantle_lithosphere < 0 || id_mid_mantle_lithosphere < 0 || 
       id_lwr_mantle_lithosphere < 0
        default_props = DensityProperties(3300.0, 3.0e-5, 1.0e-11)
        return (default_props, default_props, default_props)
    else
        return (
            DensityProperties(
                GetProperty.get_standard_density(model, id_upr_mantle_lithosphere),
                GetProperty.get_thermal_expansivity(model, id_upr_mantle_lithosphere),
                GetProperty.get_compressibility(model, id_upr_mantle_lithosphere)
            ),
            DensityProperties(
                GetProperty.get_standard_density(model, id_mid_mantle_lithosphere),
                GetProperty.get_thermal_expansivity(model, id_mid_mantle_lithosphere),
                GetProperty.get_compressibility(model, id_mid_mantle_lithosphere)
            ),
            DensityProperties(
                GetProperty.get_standard_density(model, id_lwr_mantle_lithosphere),
                GetProperty.get_thermal_expansivity(model, id_lwr_mantle_lithosphere),
                GetProperty.get_compressibility(model, id_lwr_mantle_lithosphere)
            )
        )
    end
end

function get_density_sticky_air(model::ModelData)::Float64
    matid_domains = model.materials.dicts.matid_domains
    id_air = matid_domains["Atmosphere"]
    return id_air < 0 ? 1.0 : GetProperty.get_standard_density(model, id_air)
end

function get_density_properties_continental_crust(
    model::ModelData
)::Tuple{DensityProperties,DensityProperties}
    id_upr_cont_crust, id_lwr_cont_crust = 
        DomainMaterialIDs.get_continental_crust_material_ids(model)

    if id_upr_cont_crust < 0 || id_lwr_cont_crust < 0
        default_props = DensityProperties(2800.0, 3.0e-5, 1.0e-11)
        return (default_props, default_props)
    else
        return (
            DensityProperties(
                GetProperty.get_standard_density(model, id_upr_cont_crust),
                GetProperty.get_thermal_expansivity(model, id_upr_cont_crust),
                GetProperty.get_compressibility(model, id_upr_cont_crust)
            ),
            DensityProperties(
                GetProperty.get_standard_density(model, id_lwr_cont_crust),
                GetProperty.get_thermal_expansivity(model, id_lwr_cont_crust),
                GetProperty.get_compressibility(model, id_lwr_cont_crust)
            )
        )
    end
end

function get_density_model_based_on_inputs(;
    turn_off_density_change::Bool=false,
    expansivity::Float64=3e-5,
    compressibility::Float64=1.0e-11,
    density_upper_continental_crust::Float64=2800.0,
    density_lower_continental_crust::Float64=2800.0,
    density_mantle_lithosphere::Float64=3300.0,
    density_asthenosphere::Float64=3300.0
)::DensityModel
    density_props_upr_cont_crust, density_props_lwr_cont_crust = 
        get_density_properties_continental_crust_based_on_user_input(
            density_upper_continental_crust, density_lower_continental_crust,
            expansivity, compressibility;
            turn_off_density_change=turn_off_density_change
        )

    density_props_upper_mantle_lithosphere, density_props_middle_mantle_lithosphere,
    density_props_lower_mantle_lithosphere = 
        get_density_properties_mantle_lithosphere_based_on_user_input(
            density_mantle_lithosphere, expansivity, compressibility;
            turn_off_density_change=turn_off_density_change
        )

    density_props_asthenosphere = 
        get_density_properties_asthenosphere_based_on_user_input(
            density_asthenosphere, expansivity, compressibility;
            turn_off_density_change=turn_off_density_change
        )

    return DensityModel(
        density_props_upr_cont_crust,
        density_props_lwr_cont_crust,
        density_props_upper_mantle_lithosphere,
        density_props_middle_mantle_lithosphere,
        density_props_lower_mantle_lithosphere,
        density_props_asthenosphere
    )
end

function get_density_properties_asthenosphere_based_on_user_input(
    standard_density::Float64,
    expansivity::Float64,
    compressibility::Float64;
    turn_off_density_change::Bool=true
)::DensityProperties
    if turn_off_density_change
        expansivity = 0.0
        compressibility = 0.0
    end
    return DensityProperties(standard_density, expansivity, compressibility)
end

function get_density_properties_mantle_lithosphere_based_on_user_input(
    standard_density::Float64,
    expansivity::Float64,
    compressibility::Float64;
    turn_off_density_change::Bool=true
)::Tuple{DensityProperties,DensityProperties,DensityProperties}
    if turn_off_density_change
        expansivity = 0.0
        compressibility = 0.0
    end
    props = DensityProperties(standard_density, expansivity, compressibility)
    return (props, props, props)
end

function get_density_properties_continental_crust_based_on_user_input(
    density_upper_cont_crust::Float64,
    density_lower_cont_crust::Float64,
    expansivity::Float64,
    compressibility::Float64;
    turn_off_density_change::Bool=true
)::Tuple{DensityProperties,DensityProperties}
    if turn_off_density_change
        expansivity = 0.0
        compressibility = 0.0
    end
    return (
        DensityProperties(density_upper_cont_crust, expansivity, compressibility),
        DensityProperties(density_lower_cont_crust, expansivity, compressibility)
    )
end

function get_density_props_at_y_depth(
    layer_thickness::LayerThickness,
    density_model::DensityModel,
    y_location_meters::Float64
)::Tuple{Float64,Float64,Float64}
    y_base_upr_crust, y_base_lwr_crust, y_base_upr_mantle_lithosphere,
    y_base_mid_mantle_lithosphere, y_base_lwr_mantle_lithosphere = 
        get_layer_base_depths(layer_thickness)

    if y_location_meters <= y_base_upr_crust
        return get_props(density_model.upr_cont_crust)
    elseif y_location_meters <= y_base_lwr_crust
        return get_props(density_model.lwr_cont_crust)
    elseif y_location_meters <= y_base_upr_mantle_lithosphere
        return get_props(density_model.upr_mantle_lithosphere)
    elseif y_location_meters <= y_base_mid_mantle_lithosphere
        return get_props(density_model.mid_mantle_lithosphere)
    elseif y_location_meters <= y_base_lwr_mantle_lithosphere
        return get_props(density_model.lwr_mantle_lithosphere)
    else
        return get_props(density_model.asthenosphere)
    end
end

function get_props(density_props::DensityProperties)::Tuple{Float64,Float64,Float64}
    return (
        density_props.standard_density,
        density_props.expansivity,
        density_props.compressibility
    )
end

function get_layer_base_depths(
    layer_thickness::LayerThickness
)::Tuple{Float64,Float64,Float64,Float64,Float64}
    y_base_upr_crust = layer_thickness.thick_upper_crust

    y_base_lwr_crust = layer_thickness.thick_upper_crust + layer_thickness.thick_lower_crust

    y_base_upr_mantle_lithosphere = y_base_lwr_crust + layer_thickness.thick_upper_lith

    y_base_mid_mantle_lithosphere = y_base_upr_mantle_lithosphere + 
        layer_thickness.thick_middle_lith

    y_base_lwr_mantle_lithosphere = y_base_mid_mantle_lithosphere + 
        layer_thickness.thick_lower_lith

    return (
        y_base_upr_crust,
        y_base_lwr_crust,
        y_base_upr_mantle_lithosphere,
        y_base_mid_mantle_lithosphere,
        y_base_lwr_mantle_lithosphere
    )
end

end # module