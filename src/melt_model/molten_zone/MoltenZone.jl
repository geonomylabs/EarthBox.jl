module MoltenZone

import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox: MarkerSwarm
import ..GetExtractedMoltenIDs: get_extracted_gabbroic_molten_ids
import ..GetExtractedMoltenIDs: get_mantle_molten_ids

function call_calculate_dimensions_of_molten_domain_all_drainage!(
    model::ModelData
)::Nothing
    ndrainage_basin = model.melting.parameters.extraction.ndrainage_basin.value
    xstart_drainage = model.melting.arrays.extraction.xstart_drainage.array
    xend_drainage = model.melting.arrays.extraction.xend_drainage.array

    for i in 1:ndrainage_basin
        xstart = xstart_drainage[i]
        xend = xend_drainage[i]
        (
            xmid_mol,
            ytop_mol,
            width_mol
        ) = calculate_dimensions_all_gabbroic_molten(model, xstart, xend)
        model.melting.arrays.extraction.xmid_molten_zones.array[i] = xmid_mol
        model.melting.arrays.extraction.ytop_molten_zones.array[i] = ytop_mol
        model.melting.arrays.extraction.width_molten_zones.array[i] = width_mol

        # Use partially molten mantle if no molten gabbroic material is found
        if width_mol < 0.0
            (
                xmid_mol,
                ytop_mol,
                width_mol
            ) = calculate_dimensions_mantle_molten(model, xstart, xend)
            # If no molten mantle material is found just using drainage basin
            if width_mol < 0.0
                xmid_mol = (xend - xstart)/2.0
                ytop_mol = 10_000.0
                width_mol = xend - xstart
            end
        end

        print_info = false
        if print_info
            print_molten_zone_info(i, xmid_mol, ytop_mol, width_mol)
        end
    end
    return nothing
end

function print_molten_zone_info(
    i::Int,
    xmid_mol::Float64,
    ytop_mol::Float64,
    width_mol::Float64
)::Nothing
    print_info(">> Molten zone dimensions for drainage basin $(i)", level=1)
    print_info("xmid_mol: $(xmid_mol)", level=2)
    print_info("ytop_mol: $(ytop_mol)", level=2)
    print_info("width_mol: $(width_mol)", level=2)
    print_info("", level=1)
    return nothing
end

""" Calculate dimensions of molten zone.

# Arguments
- model::ModelData
    - Model data structure.
- xstart::Float64
    - Starting x-coordinate of drainage basin.
- xend::Float64
    - Ending x-coordinate of drainage basin.

# Returns
- xmid_molten::Float64
    - X-location of the middle of molten zone in meters.
- ytop_molten::Float64
    - Minimum y-location of the molten zone in meters.
- width_molten::Float64
    - Width of molten zone in meters.
"""
function calculate_dimensions_all_gabbroic_molten(
    model::ModelData,
    xstart::Float64 = -1e39,
    xend::Float64 = 1e39
)::Tuple{Float64, Float64, Float64}
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    markers_matid = model.markers.arrays.material.marker_matid.array
    extracted_molten_matids = get_extracted_gabbroic_molten_ids(model)

    (
        xmid_molten, ytop_molten, width_molten
    ) = MarkerSwarm.calculate_dimensions_general(
        marker_x, marker_y, markers_matid,
        extracted_molten_matids,
        xstart, xend
    )

    return xmid_molten, ytop_molten, width_molten
end

""" Calculate dimensions of mantle molten zone.

# Arguments
- model::ModelData
    - Model data structure.
- xstart::Float64
    - Starting x-coordinate of drainage basin.
- xend::Float64
    - Ending x-coordinate of drainage basin.

# Returns
- xmid_molten::Float64
    - X-location of the middle of molten zone in meters.
- ytop_molten::Float64
    - Minimum y-location of the molten zone in meters.
- width_molten::Float64
    - Width of molten zone in meters.
"""
function calculate_dimensions_mantle_molten(
    model::ModelData,
    xstart::Float64 = -1e39,
    xend::Float64 = 1e39
)::Tuple{Float64, Float64, Float64}
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    markers_matid = model.markers.arrays.material.marker_matid.array
    extracted_molten_matids = get_mantle_molten_ids(model)

    (
        xmid_molten, ytop_molten, width_molten
    ) = MarkerSwarm.calculate_dimensions_general(
        marker_x, marker_y, markers_matid,
        extracted_molten_matids,
        xstart, xend
    )

    return xmid_molten, ytop_molten, width_molten
end

function get_x_limits_of_molten_zone(
    model::ModelData,
    idrainage_basin::Int64
)::Tuple{Float64, Float64}
    width_molten_zones = model.melting.arrays.extraction.width_molten_zones.array
    xmid_molten_zones = model.melting.arrays.extraction.xmid_molten_zones.array

    width_molten_domain = width_molten_zones[idrainage_basin]
    xmid_molten_domain = xmid_molten_zones[idrainage_basin]
    if width_molten_domain > 0.0
        xmin_molten_domain = xmid_molten_domain - width_molten_domain/2
        xmax_molten_domain = xmid_molten_domain + width_molten_domain/2
    else
        xmin_molten_domain = -1e39
        xmax_molten_domain = -1e39
    end
    return xmin_molten_domain, xmax_molten_domain
end
    
end # module 