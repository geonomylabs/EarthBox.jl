module MarkerCompaction

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_swarm
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_swarm_opt
import EarthBox: MathTools
import ..CompactionTools: compact_or_decompact

""" Compact sediment and advect markers given deposited sediment thickness.

The initial sediment thickness may include thickness not yet populated with 
sediment markers. For these cases default porosity and decay depth values are 
used to calculate compaction properties.

Sticky markers are advected using the compaction displacement at the top of 
the sediment markers and a linear change in displacement from the 
sediment-sticky interface to the top of the sticky layer.

Sticky markers within the domain of newly deposited sediment are converted to 
sediment by a separate function.

Updated Arrays
-------------
markers_y : Vector{Float64}
    The y-coordinates of the markers.

markers_max_burial_depth : Vector{Float64}
    The current maximum burial depth of the markers (meters).

Algorithm
---------
1. Make the compaction array.
2. Calculate vertical cell geometry for compaction meshes.
3. Assign the topo x index to the markers.
4. Assign the compaction y index to the markers and calculate unit distance
    from cell top.
5. Calculate the compaction properties for the compaction array using
    marker properties.
6. Compact sediment mesh.
7. Calculate total mesh compaction y-displacement.
8. Calculate total marker compaction displacement using mesh displacement.
9. Apply compaction displacement to sediment markers.
10. Apply compaction displacement to sticky markers
"""
function compact_sediment_and_advect_markers(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    sediment_thickness_initial_gridx::Vector{Float64},
    sticky_thickness_initial_gridx::Vector{Float64},
    new_sediment_thickness_gridx::Vector{Float64},
    marker_porosity_initial::Vector{Float32},
    marker_decay_depth::Vector{Float32},
    markers_max_burial_depth::Vector{Float32},
    markers_indices_sedimentary_basin::Vector{Int64},
    markers_indices_sticky::Vector{Int64},
    x_sorted_marker_indices_sedimentary_basin::Vector{Int64},
    markers_x::Vector{Float64},
    markers_y::Vector{Float64},
    search_radius::Float64,
    nsmooth_top_bottom::Int;
    default_porosity_initial::Float64 = 0.4,
    default_decay_depth::Float64 = 2000.0
)::Vector{Float64}
    compaction_array = make_compaction_array(topo_gridx)
    calculate_compaction_meshes(
        topo_gridx, topo_gridy,
        sediment_thickness_initial_gridx, compaction_array
    )
    markers_topo_xindex = assign_topo_xindex_to_markers(
        topo_gridx, markers_x, markers_indices_sedimentary_basin
    )
    (
        markers_compaction_yindex,
        markers_unit_distance_from_cell_top
    ) = assign_compaction_yindex_to_markers(
        markers_y, markers_indices_sedimentary_basin,
        compaction_array, markers_topo_xindex
    )
    calculate_compaction_properties_for_compaction_array_opt(
        compaction_array, markers_topo_xindex, markers_compaction_yindex,
        marker_porosity_initial, marker_decay_depth, markers_max_burial_depth,
        markers_indices_sedimentary_basin
    )
    calculate_average_compaction_properties_in_compaction_cells(
        compaction_array;
        default_porosity_initial=default_porosity_initial,
        default_decay_depth=default_decay_depth
    )
    compact_sediment_mesh(
        new_sediment_thickness_gridx, topo_gridy, compaction_array
    )
    sediment_thickness_initial_compacted_gridx = 
        calculate_mesh_compacted_thickness(compaction_array)
    calculate_total_mesh_compaction_displacement(compaction_array)
    calculate_updated_mesh_burial_depth(compaction_array)
    update_marker_max_burial(
        markers_max_burial_depth, markers_topo_xindex,
        markers_compaction_yindex, markers_unit_distance_from_cell_top,
        markers_indices_sedimentary_basin,
        compaction_array
    )
    total_marker_compaction_displacement = 
        calculate_total_marker_compaction_displacement(
            markers_y, markers_topo_xindex, markers_compaction_yindex,
            markers_unit_distance_from_cell_top,
            markers_indices_sedimentary_basin,
            compaction_array
        )
    (
        top_marker_sediment_pre_compaction, _
    ) = calculate_top_and_bottom_of_swarm_opt(
        x_sorted_marker_indices_sedimentary_basin,
        markers_x, markers_y, topo_gridx,
        search_radius;
        use_smoothing=true,
        nsmooth=nsmooth_top_bottom
    )
    apply_compaction_displacement_to_sediment_markers(
        markers_indices_sedimentary_basin,
        markers_y,
        total_marker_compaction_displacement,
        markers_compaction_yindex
    )
    (
        top_marker_sediment_post_compaction, _
    ) = calculate_top_and_bottom_of_swarm_opt(
        x_sorted_marker_indices_sedimentary_basin,
        markers_x, markers_y, topo_gridx,
        search_radius;
        use_smoothing=true,
        nsmooth=nsmooth_top_bottom
    )
    max_sticky_displacement = (
        top_marker_sediment_post_compaction -
        top_marker_sediment_pre_compaction
    )
    apply_compaction_displacement_to_sticky_markers(
        markers_x,
        markers_y,
        topo_gridx,
        markers_indices_sticky,
        sticky_thickness_initial_gridx,
        max_sticky_displacement
    )
    return sediment_thickness_initial_compacted_gridx
end

function print_compaction_info(
    sediment_thickness_initial_gridx::Vector{Float64},
    sediment_thickness_initial_compacted_gridx::Vector{Float64},
    total_marker_compaction_displacement::Vector{Float64},
    max_sticky_displacement::Vector{Float64}
)::Nothing
    """ Print compaction information.
    """
    println(
        "Max of initial sediment thickness (m) ",
        maximum(sediment_thickness_initial_gridx)
    )

    println(
        "Max of compacted sediment thickness (m) ",
        maximum(sediment_thickness_initial_compacted_gridx)
    )

    println(
        "Max of total marker compaction displacement (m) ",
        maximum(total_marker_compaction_displacement)
    )
    println(
        "Max of total sticky displacement (m) ",
        maximum(max_sticky_displacement)
    )
end

function calculate_swarm_indices_for_sediment_and_sticky(
    model::ModelData
)::Tuple{Vector{Int64}, Vector{Int64}}
    markers_x = model.markers.arrays.location.marker_x.array
    markers_matid = model.markers.arrays.material.marker_matid.array

    sedimentary_basin_ids = get_sedimentary_basin_material_ids(
        model.materials.dicts.matid_types
    )
    markers_indices_sedimentary_basin = calculate_marker_swarm_indices(
        markers_x, markers_matid, sedimentary_basin_ids
    )

    sticky_ids = get_sticky_material_ids(
        model.materials.dicts.matid_types
    )
    markers_indices_sticky = calculate_marker_swarm_indices(
        markers_x, markers_matid, sticky_ids
    )

    return markers_indices_sedimentary_basin, markers_indices_sticky
end

function calculate_x_sorted_swarm_indices(
    model::ModelData,
    marker_indices_swarm::Vector{Int64}
)::Vector{Int64}
    marker_x = model.markers.arrays.location.marker_x.array
    # Sort markers by x-coordinate
    sorted_indices = sortperm(marker_x[marker_indices_swarm])
    sorted_marker_indices_swarm = marker_indices_swarm[sorted_indices]
    return sorted_marker_indices_swarm
end

function calculate_x_sorted_swarm_indices_from_marker_x(
    marker_x::Vector{Float64},
    marker_indices_swarm::Vector{Int64}
)::Vector{Int64}
    # Sort markers by x-coordinate
    sorted_indices = sortperm(marker_x[marker_indices_swarm])
    sorted_marker_indices_swarm = marker_indices_swarm[sorted_indices]
    return sorted_marker_indices_swarm
end

function get_sticky_material_ids(
    types::Dict{String, Vector{Int16}}
)::Vector{Int16}
    if length(types["StickyAir"]) > 0
        matid_sticky_air = types["StickyAir"][1]
    else
        matid_sticky_air = Int16(-1)
    end
    if length(types["StickyWater"]) > 0
        matid_sticky_water = types["StickyWater"][1]
    else
        matid_sticky_water = Int16(-1)
    end
    sticky_ids = [matid_sticky_air, matid_sticky_water]
    return sticky_ids
end

function get_sedimentary_basin_material_ids(
    types::Dict{String, Vector{Int16}}
)::Vector{Int16}
    if length(types["Sediment"]) > 0
        matid_sediment = types["Sediment"][1]
    else
        matid_sediment = Int16(-1)
    end
    if length(types["SolidifiedBasalt"]) > 0
        matid_solidified_basalt = types["SolidifiedBasalt"][1]
    else
        matid_solidified_basalt = Int16(-1)
    end
    if length(types["Salt"]) > 0
        matid_salt = types["Salt"][1]
    else
        matid_salt = Int16(-1)
    end
    sedimentary_basin_ids = [matid_sediment, matid_solidified_basalt, matid_salt]
    return sedimentary_basin_ids
end

function calculate_marker_swarm_indices(
    marker_x::Vector{Float64},
    marker_matid::Vector{Int16},
    swarm_material_ids::Vector{Int16}
)::Vector{Int64}
    """ Calculate marker indices for a swarm with a set of id's.

    Inputs
    ------
    marker_x: Vector{Float64}
        The x-coordinates of the markers.

    marker_matid: Vector{Int16}
        The material ids of the markers.

    swarm_material_ids: Vector{Int16}
        The material ids of the swarm.

    Returns
    -------
    marker_swarm_indices : Vector{Int64}
        The index of the sedimentary basin markers in the original marker arrays.
    """
    nmarkers = length(marker_x)
    temp_marker_index = ones(Int64, nmarkers)

    nmarkers_swarm = 0
    for i in 1:nmarkers
        if marker_matid[i] in swarm_material_ids
            temp_marker_index[nmarkers_swarm+1] = i
            nmarkers_swarm += 1
        end
    end

    marker_swarm_indices = temp_marker_index[1:nmarkers_swarm]
    return marker_swarm_indices
end

function make_compaction_array(
    topo_gridx::Vector{Float64}
)::Array{Float64, 3}
    """ Make the compaction array.

    Returns
    -------
    compaction_array : Array{Float64, 3}
        The compaction array:
        - The first index is the x index of the compaction grid equivalent to the 
          topography grid.
        - The second index is the cell index in the vertical y-direction of the 
          compaction grid.
        - The third index refers to the following properties:
            - 1: y coordinate of the compaction grid.
            - 2: Initial porosity in fraction (from markers).
            - 3: Porosity decay depth in meters (from markers).
            - 4: Number of markers in cell.
            - 5: Thickness of layer in meters.
            - 6: Incremental thickness change of cell (meters).
            - 7: Total compaction y-displacement relative to top of un-compacted 
                 sediment (meters).
            - 8: Maximum burial depth of cell (meters) (from markers).
            - 9: Updated burial depth of cell (meters).
    """
    ny_comp = 20
    compaction_array = zeros(Float64, length(topo_gridx), ny_comp, 9)
    return compaction_array
end

""" Calculate vertical cell geometry for compaction meshes.

Vertical compaction meshes are hung on the topography grid and have a
thickness equal to sediment thickness.

Updated Arrays
-------------
compaction_array : Array{Float64, 3}
    The compaction array. Updated components are as follows:
    - 1: y coordinate of the tops pf compaction grid cells.
    - 5: Thickness of layer in meters.
"""
function calculate_compaction_meshes(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    sediment_thickness_gridx::Vector{Float64},
    compaction_array::Array{Float64, 3}
)::Nothing
    ncells = size(compaction_array, 2)
    ntopo = length(topo_gridx)
    for j in 1:ncells
        for i in 1:ntopo
            y_topo = topo_gridy[i]
            if sediment_thickness_gridx[i] > 0.0
                dy = sediment_thickness_gridx[i] / Float64(ncells)
            else
                dy = 0.0
            end
            compaction_array[i, j, 1] = y_topo + Float64(j-1) * dy # depth y-coordinate
            compaction_array[i, j, 5] = dy # thickness of cell
        end
    end
end

function assign_topo_xindex_to_markers(
    topo_gridx::Vector{Float64},
    markers_x::Vector{Float64},
    marker_indices_sedimentary_basin::Vector{Int64}
)::Vector{Int64}
    nswarm = length(marker_indices_sedimentary_basin)
    nmarkers = length(markers_x)
    ntopo = length(topo_gridx)
    markers_topo_xindex = Vector{Int64}(undef, nmarkers)
    Threads.@threads for i in 1:nswarm
        imarker = marker_indices_sedimentary_basin[i]
        x_marker = markers_x[imarker]
        xindex = 1 # the original code had xindex = 0
        for j in 1:ntopo
            if x_marker <= topo_gridx[j]
                xindex = j
                break
            end
        end
        markers_topo_xindex[imarker] = xindex
    end
    return markers_topo_xindex
end

""" Assign the compaction y index to the markers and calculate unit distance
    from cell top.

Returns
-------
markers_compaction_yindex : Vector{Int64}
    The y index of the compaction grid cell. -1 indicates that the marker
    is not in the compaction grid.

markers_unit_distance_from_cell_top : Vector{Float64}
    The unit distance from the top of the compaction grid cell.
"""
function assign_compaction_yindex_to_markers(
    markers_y::Vector{Float64},
    marker_indices_sedimentary_basin::Vector{Int64},
    compaction_array::Array{Float64, 3},
    markers_topo_xindex::Vector{Int64}
)::Tuple{Vector{Int64}, Vector{Float64}}
    markers_compaction_yindex = Vector{Int64}(undef, length(markers_y))
    markers_unit_distance_from_cell_top = Vector{Float64}(undef, length(markers_y))
    ncells = size(compaction_array, 2)
    nmarkers_sediment = length(marker_indices_sedimentary_basin)
    Threads.@threads for i in 1:nmarkers_sediment
        imarker = marker_indices_sedimentary_basin[i]
        xindex = markers_topo_xindex[imarker]
        y_marker = markers_y[imarker]
        compaction_yindex = -1
        unit_distance = 0.0
        for j in 1:ncells
            dy = compaction_array[xindex, j, 5]
            y_top_cell = compaction_array[xindex, j, 1]
            y_bottom_cell = y_top_cell + dy
            if y_top_cell <= y_marker < y_bottom_cell
                compaction_yindex = j
                if dy > 0.0
                    unit_distance = 1.0 - (y_marker - y_top_cell) / dy
                else
                    unit_distance = 0.0
                end
                break
            end
        end
        markers_compaction_yindex[imarker] = compaction_yindex
        markers_unit_distance_from_cell_top[imarker] = unit_distance
    end
    return markers_compaction_yindex, markers_unit_distance_from_cell_top
end

""" Calculate the compaction properties for the compaction array.

Updated Arrays
-------------
compaction_array : Array{Float64, 3}
    The compaction array. The following components are updated:
    - 2: Initial porosity in fraction.
    - 3: Porosity decay depth in meters.
    - 4: Number of markers in cell (used to calculate averages).
    - 8: Maximum burial depth of cell (meters).
"""
function calculate_compaction_properties_for_compaction_array_opt(
    compaction_array::Array{Float64, 3},
    markers_topo_xindex::Vector{Int64},
    markers_compaction_yindex::Vector{Int64},
    marker_porosity_initial::Vector{Float32},
    marker_decay_depth::Vector{Float32},
    markers_max_burial_depth::Vector{Float32},
    marker_indices_sedimentary_basin::Vector{Int64}
)::Nothing
    ntopo = size(compaction_array, 1)
    ncells = size(compaction_array, 2)
    nmarkers_sediment = length(marker_indices_sedimentary_basin)

    Threads.@threads for k in 1:nmarkers_sediment
        imarker = marker_indices_sedimentary_basin[k]
        compaction_yindex = markers_compaction_yindex[imarker]
        if !(compaction_yindex < 1)
            topo_xindex = markers_topo_xindex[imarker]
            if 1 <= topo_xindex <= ntopo && 1 <= compaction_yindex <= ncells
                compaction_array[topo_xindex, compaction_yindex, 2] += 
                    marker_porosity_initial[imarker]
                compaction_array[topo_xindex, compaction_yindex, 3] += 
                    marker_decay_depth[imarker]
                compaction_array[topo_xindex, compaction_yindex, 8] = 
                    markers_max_burial_depth[imarker]
                compaction_array[topo_xindex, compaction_yindex, 4] += 1.0
            end
        end
    end
end

""" Finalize the compaction properties for the compaction array.

Updated Arrays
-------------
compaction_array : Array{Float64, 3}
    The compaction array. The following components are updated:
    - 2: Initial porosity in fraction.
    - 3: Porosity decay depth in meters.
    - 8: Maximum burial depth of cell (meters).
"""
function calculate_average_compaction_properties_in_compaction_cells(
    compaction_array::Array{Float64, 3};
    default_porosity_initial::Float64 = 0.4,
    default_decay_depth::Float64 = 2000.0
)::Nothing
    ntopo = size(compaction_array, 1)
    ncells = size(compaction_array, 2)
    
    Threads.@threads for j in 1:ncells
        for i in 1:ntopo
            if compaction_array[i, j, 4] > 0.0
                compaction_array[i, j, 2] /= compaction_array[i, j, 4]
                compaction_array[i, j, 3] /= compaction_array[i, j, 4]
                compaction_array[i, j, 8] /= compaction_array[i, j, 4]
            else
                compaction_array[i, j, 2] = default_porosity_initial
                compaction_array[i, j, 3] = default_decay_depth
                compaction_array[i, j, 8] = 0.0
            end
        end
    end
end

""" Compact sediment mesh.

For each compaction mesh associated with x-index i, calculate new
thicknesses of layers by looping through each vertical cell j starting from
the top. Cell thickness is compacted only if the new submud depth of
the cell top is greater than the maximum burial depth of the cell based
on average maximum burial depths from all markers within cell.

Updated Arrays
-------------
compaction_array : Array{Float64, 3}
    The compaction array. Updated components are as follows:
    - 1: y coordinate of the tops pf compaction grid cells.
    - 5: Thickness of layer in meters.
    - 6: Incremental thickness change of cell (meters).
"""
function compact_sediment_mesh(
    new_sediment_thickness_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    compaction_array::Array{Float64, 3}
)::Nothing
    ntopo = size(compaction_array, 1)
    ncells = size(compaction_array, 2)
    
    # The structure of the array may need to be modified to optimize memory access
    # for 1-based arrays
    Threads.@threads for i in 1:ntopo
        new_sediment_thickness = new_sediment_thickness_gridx[i]
        if new_sediment_thickness > 0.0
            for j in 1:ncells
                top_initial_submud = compaction_array[i, j, 1] - topo_gridy[i]
                thickness_initial = compaction_array[i, j, 5]
                bottom_initial = top_initial_submud + thickness_initial

                if j == 1
                    top_new_submud = new_sediment_thickness
                else
                    top_above_submud = compaction_array[i, j-1, 1] - topo_gridy[i]
                    thickness_above = compaction_array[i, j-1, 5]
                    top_new_submud = top_above_submud + thickness_above
                end

                porosity_initial = compaction_array[i, j, 2]
                depth_decay_term = 1.0 / compaction_array[i, j, 3]
                average_marker_max_burial_depth = compaction_array[i, j, 8]

                # Only compact if base of cell is deeper than average max
                # burial depth of markers in cell
                if bottom_initial > average_marker_max_burial_depth
                    thickness_new = compact_or_decompact(
                        porosity_initial, depth_decay_term,
                        top_initial_submud, bottom_initial, top_new_submud
                    )
                else
                    thickness_new = thickness_initial
                end

                # Update total y-depth of top of cell
                compaction_array[i, j, 1] = top_new_submud + topo_gridy[i]
                # Update thickness of cell
                compaction_array[i, j, 5] = thickness_new
                # Update incremental thickness change of cell
                compaction_array[i, j, 6] = thickness_new - thickness_initial
            end
        end
    end
end

""" Calculate total compacted sediment thickness using mesh.

Returns
-------
sediment_thickness_compacted : Vector{Float64}
    The total compacted sediment thickness using mesh.
"""
function calculate_mesh_compacted_thickness(
    compaction_array::Array{Float64, 3}
)::Vector{Float64}
    ntopo = size(compaction_array, 1)
    ncells = size(compaction_array, 2)
    sediment_thickness_compacted = zeros(Float64, ntopo)
    
    # Optimization needed for 1-based arrays
    Threads.@threads for i in 1:ntopo
        for j in 1:ncells
            sediment_thickness_compacted[i] += abs(compaction_array[i, j, 5])
        end
    end
    return sediment_thickness_compacted
end

""" Calculate total compaction y-displacement

Total displacement is the sum of incremental thickness changes of the cell
and all cells below the cell. Displacement is calculated relative to the
basement.

Updated Arrays
-------------
compaction_array : Array{Float64, 3}
    The compaction array. Updated components are as follows:
    - 7: Total compaction y-displacement relative to top of decompacted
            sediment (meters).
    The following component is used in the calculation:
    - 6: Incremental thickness change of cell (meters).
"""
function calculate_total_mesh_compaction_displacement(
    compaction_array::Array{Float64, 3}
)::Nothing
    ntopo = size(compaction_array, 1)
    ncells = size(compaction_array, 2)
    Threads.@threads for i in 1:ntopo
        for j in 1:ncells
            for k in j:ncells
                compaction_array[i, j, 7] += abs(compaction_array[i, k, 6])
            end
        end
    end
end

""" Calculate the maximum burial depth of the mesh.

Burial depth is the sum of all cell thicknesses above a given cell.

Updated Arrays
-------------
compaction_array : Array{Float64, 3}
    The compaction array. Updated components are as follows:
    - 9: Updated burial depth of cell (meters).
    The following component is used in the calculation:
    - 5: Thickness of layer in meters.
"""
function calculate_updated_mesh_burial_depth(
    compaction_array::Array{Float64, 3}
)::Nothing
    ntopo = size(compaction_array, 1)
    ncells = size(compaction_array, 2)
    
    # Optimization needed for 1-based arrays
    Threads.@threads for i in 1:ntopo
        for j in 1:ncells
            for k in 1:j-1
                compaction_array[i, j, 9] += compaction_array[i, k, 5]
            end
        end
    end
end

""" Update marker maximum burial depth.

This function updates the maximum burial depth of the markers using the
the current burial depth of the compaction mesh.

Updated Arrays
-------------
markers_max_burial_depth : Vector{Float64}
    The current maximum burial depth of the markers (meters).
"""
function update_marker_max_burial(
    markers_max_burial_depth::Vector{Float32},
    markers_topo_xindex::Vector{Int64},
    markers_compaction_yindex::Vector{Int64},
    markers_unit_distance_from_cell_top::Vector{Float64},
    markers_indices_sedimentary_basin::Vector{Int64},
    compaction_array::Array{Float64, 3}
)::Nothing
    nswarm = length(markers_indices_sedimentary_basin)
    ncells = size(compaction_array, 2)
    
    for i in 1:nswarm
        imarker = markers_indices_sedimentary_basin[i]
        xindex = markers_topo_xindex[imarker]
        cell_index = markers_compaction_yindex[imarker]
        if cell_index != -1
            unit_distance = markers_unit_distance_from_cell_top[imarker]
            depth_top = compaction_array[xindex, cell_index, 9]
            if cell_index == ncells
                depth_bottom = depth_top + 1.0
            else
                depth_bottom = compaction_array[xindex, cell_index + 1, 9]
            end
            current_marker_burial_depth = (
                depth_top +
                unit_distance * (depth_bottom - depth_top)
            )
            markers_max_burial_depth[imarker] = max(
                markers_max_burial_depth[imarker],
                current_marker_burial_depth
            )
        end
    end
end

""" Calculate the total marker compaction displacement.

Returns
-------
total_marker_compaction_displacement : Vector{Float64}
    The total marker compaction displacement.
"""
function calculate_total_marker_compaction_displacement(
    markers_y::Vector{Float64},
    markers_topo_xindex::Vector{Int64},
    markers_compaction_yindex::Vector{Int64},
    markers_unit_distance_from_cell_top::Vector{Float64},
    markers_indices_sedimentary_basin::Vector{Int64},
    compaction_array::Array{Float64, 3}
)::Vector{Float64}
    nswarm = length(markers_indices_sedimentary_basin)
    nmarkers = length(markers_y)
    total_marker_compaction_displacement = zeros(Float64, nmarkers)
    ncells = size(compaction_array, 2)
    
    Threads.@threads for i in 1:nswarm
        imarker = markers_indices_sedimentary_basin[i]
        xindex = markers_topo_xindex[imarker]
        cell_index = markers_compaction_yindex[imarker]
        if !(cell_index < 1)
            unit_distance = markers_unit_distance_from_cell_top[imarker]
            displacement_top = compaction_array[xindex, cell_index, 7]
            if cell_index == ncells
                displacement_bottom = 0.0
            else
                displacement_bottom = compaction_array[xindex, cell_index + 1, 7]
            end
            total_marker_compaction_displacement[imarker] = (
                displacement_top +
                unit_distance * (displacement_bottom - displacement_top)
            )
        end
    end
    return total_marker_compaction_displacement
end

""" Apply compaction displacement to sediment markers.

Updated Arrays
-------------
markers_y : Vector{Float64}
    The y-coordinates of the markers.
"""
function apply_compaction_displacement_to_sediment_markers(
    marker_indices_sedimentary_basin::Vector{Int64},
    markers_y::Vector{Float64},
    total_marker_compaction_displacement::Vector{Float64},
    markers_compaction_yindex::Vector{Int64}
)::Nothing
    nmarkers_sediment = length(marker_indices_sedimentary_basin)
    for i in 1:nmarkers_sediment
        imarker = marker_indices_sedimentary_basin[i]
        yindex = markers_compaction_yindex[imarker]
        y_initial = markers_y[imarker]
        if !(yindex < 1)
            markers_y[imarker] = (
                y_initial + total_marker_compaction_displacement[imarker]
            )
        end
    end
end

""" Apply compaction displacement to sticky markers.

Algorithm
---------
1. Calculate compaction displacement factors for sticky markers.
2. Calculate sticky marker displacement.
3. Move sticky markers using compaction field.
"""
function apply_compaction_displacement_to_sticky_markers(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    topo_gridx::Vector{Float64},
    marker_indices_sticky::Vector{Int64},
    sticky_thickness::Vector{Float64},
    max_sticky_displacement::Vector{Float64}
)::Nothing
    marker_displacement_factors = calculate_sticky_compaction_displacement_factors_opt(
        marker_x, marker_y, topo_gridx,
        marker_indices_sticky, sticky_thickness
    )
    marker_displacement = calculate_sticky_marker_displacement_opt(
        marker_x, topo_gridx, marker_displacement_factors,
        max_sticky_displacement, marker_indices_sticky
    )
    move_sticky_markers_using_compaction_field(
        marker_y, marker_displacement, marker_indices_sticky
    )
end

""" Calculate the maximum compaction displacement field at each x-grid cell 
along the sediment-sticky interface.

Returns
-------
compaction_displacement_max : Vector{Float64}
    The maximum compaction displacement field at each x-grid cell along the 
    sediment-sticky interface.
"""
function calculate_compaction_displacement_max(
    compaction_array::Array{Float64, 3}
)::Vector{Float64}
    ntopo = size(compaction_array, 1)
    compaction_displacement_max = zeros(Float64, ntopo)
    for i in 1:ntopo
        compaction_displacement_max[i] = compaction_array[i, 1, 7]
    end
    return compaction_displacement_max
end

""" Calculate compaction displacement factors for sticky air/water.

Displacement factors are unit vectors that indicate the displacement of
sticky markers relative to the sediment water interface due to
compaction. Vertical displacement is maximum at the interface between
sticky air/water and sediment.

Inputs
------
marker_x: Vector{Float64}
    The x-coordinates of the markers.

marker_y: Vector{Float64}
    The y-coordinates of the markers.

topo_gridx: Vector{Float64}
    The x-coordinates of the topography grid.

marker_indices_sticky: Vector{Int64}
    Indices of sticky markers.

sticky_thickness_markers: Vector{Float64}
    Initial sticky thickness markers at each x-grid cell.

Returns
-------
marker_displacement_factors: Vector{Float64}
    The compaction displacement factors for sticky markers.
"""
function calculate_sticky_compaction_displacement_factors_opt(
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    topo_gridx::Vector{Float64},
    marker_indices_sticky::Vector{Int64},
    sticky_thickness_markers::Vector{Float64}
)::Vector{Float64}
    nswarm = length(marker_indices_sticky)
    marknum = length(marker_x)
    marker_displacement_factors = Vector{Float64}(undef, marknum)
    Threads.@threads for i in 1:nswarm
        imarker = marker_indices_sticky[i]
        x_marker = marker_x[imarker]
        y_marker = marker_y[imarker]
        idx = searchsortedfirst(topo_gridx, x_marker) - 1
        if idx < 1
            idx = 1
        elseif idx >= length(topo_gridx) - 1
            idx = length(topo_gridx) - 2
        end
        x0, x1 = topo_gridx[idx], topo_gridx[idx + 1]
        y0, y1 = sticky_thickness_markers[idx], sticky_thickness_markers[idx + 1]
        if x1 != x0
            sticky_thickness = y0 + (y1 - y0) * (x_marker - x0) / (x1 - x0)
        else
            sticky_thickness = y0
        end
        if abs(sticky_thickness) != 0.0
            displacement_factor = y_marker / sticky_thickness
        else
            displacement_factor = 0.0
        end
        marker_displacement_factors[imarker] = displacement_factor
    end
    return marker_displacement_factors
end

""" Calculate marker displacement.

Inputs
------
markers_x: Vector{Float64}
    The x-coordinates of the markers.

topo_gridx: Vector{Float64}
    The x-coordinates of the topography grid.

marker_displacement_factors: Vector{Float64}
    Displacement factors for sticky markers.

compaction_displacement_max: Vector{Float64}
    Maximum compaction displacement field at each x-grid cell along the
    sediment-sticky interface.

Returns
-------
marker_displacement: Vector{Float64}
    The marker displacement.
"""
function calculate_sticky_marker_displacement_opt(
    markers_x::Vector{Float64},
    topo_gridx::Vector{Float64},
    marker_displacement_factors::Vector{Float64},
    compaction_displacement_max::Vector{Float64},
    marker_indices_sticky::Vector{Int64}
)::Vector{Float64}
    marknum = length(markers_x)
    nswarm = length(marker_indices_sticky)
    marker_displacement = Vector{Float64}(undef, marknum)
    Threads.@threads for i in 1:nswarm
        imarker = marker_indices_sticky[i]
        x_marker = markers_x[imarker]
        displacement_factor = marker_displacement_factors[imarker]
        displacement_max = MathTools.linear_interp_bisection(
            topo_gridx, compaction_displacement_max, x_marker
        )
        displacement = displacement_factor * displacement_max
        marker_displacement[imarker] = displacement
    end

    return marker_displacement
end

""" Calculate marker displacement.

Updated Arrays
-------------
markers_y: Vector{Float64}
    The y-coordinates of the markers.

Inputs
------
markers_y: Vector{Float64}
    The y-coordinates of the markers.

marker_displacement: Vector{Float64}
    Displacement of sticky markers due to compaction.
"""
function move_sticky_markers_using_compaction_field(
    markers_y::Vector{Float64},
    marker_displacement::Vector{Float64},
    marker_indices_sticky::Vector{Int64}
)::Nothing
    nswarm = length(marker_indices_sticky)
    Threads.@threads for i in 1:nswarm
        imarker = marker_indices_sticky[i]
        markers_y[imarker] = (
            markers_y[imarker] + marker_displacement[imarker]
        )
    end
end

end # module 