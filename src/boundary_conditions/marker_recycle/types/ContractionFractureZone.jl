module ContractionFractureZone

import EarthBox.ModelDataContainer: ModelData
import ..ContractionRecycleLocations: calculate_recycled_locations
import ..Coordinates: reset_recycled_marker_coordinates
import ..ResetMarkers: reset_recycled_marker_properties!
import ..InitializeMarkerRecycling: initialize_recycling
import ..FractureZoneProperties: get_subsurface_reset_properties_fracture_zone
import ..CheckDicts: check_domain_dict

function recycle_markers_symmetric!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    check_required_types_and_domains_fracture_zone(model)
    execute_recycle_steps(model, inside_flags; use_asymmetric_contraction=false)
    return nothing
end

function recycle_markers_asymmetric!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    check_required_types_and_domains_fracture_zone(model)
    execute_recycle_steps(model, inside_flags; use_asymmetric_contraction=true)
    return nothing
end

function check_required_types_and_domains_fracture_zone(model::ModelData)::Nothing
    matid_types = model.materials.dicts.matid_types
    type_dict = Dict(
        "StickyAir" => matid_types["StickyAir"][1],
        "StickyWater" => matid_types["StickyWater"][1]
    )
    check_domain_dict(type_dict, "type", "contractional")

    domains = model.materials.dicts.matid_domains
    domain_dict = Dict(
        "Asthenosphere" => domains["Asthenosphere"],
        "SedimentaryBasin" => domains["SedimentaryBasin"],
        "BasalticOceanicCrust" => domains["BasalticOceanicCrust"],
        "GabbroicOceanicCrust" => domains["GabbroicOceanicCrust"],
        "MantleLithosphere" => domains["MantleLithosphere"],
        "WeakMantleLithosphere" => domains["WeakMantleLithosphere"]
    )
    check_domain_dict(domain_dict, "domain", "contractional")
    return nothing
end

function execute_recycle_steps(
    model::ModelData, 
    inside_flags::Vector{Int8};
    use_asymmetric_contraction::Bool = false
)::Nothing
    (
        outside_indices, nrecycle_air, nrecycle_subsurface
    ) = initialize_recycling(model, inside_flags)
    recycled_coordinates = calculate_recycled_locations(
        model, nrecycle_air, nrecycle_subsurface, use_asymmetric_contraction)
    marker_recycling = reset_recycled_marker_coordinates(
        model, recycled_coordinates, outside_indices)
    subsurface_props = get_subsurface_reset_properties_fracture_zone(
        model, marker_recycling)
    reset_recycled_marker_properties!(model, marker_recycling, subsurface_props)
    return nothing
end

end # module 