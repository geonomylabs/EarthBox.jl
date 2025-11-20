module MarkerAdvection

import Plots
import EarthBox.ModelDataContainer: ModelData
import EarthBox.SedimentWaterInterface: get_depth
import EarthBox.ModelStructureManager.TopAndBottom: calculate_layer_thickness
import EarthBox.MathTools: linear_interp_at_x_location
import EarthBox.ModelStructureManager.TopAndBottom: calculate_layer_thickness
import EarthBox.DataStructures: SedimentTransportParameters

""" Advect markers using compaction displacement field.

This function must be called prior to updating the y-coordinate of
topography stored in gridt[i, 1] since the pre-transport mudline is used
as a reference to calculate the compaction displacement field.

# Updated Arrays
- `marker_y`: Array of marker y-coordinates (meters)
- `markers_max_burial_depth`: Array of marker maximum burial depths (meters)

# Arguments
- `model`: Model data
- `sediment_thickness_markers_initial`: Initial sediment thickness markers at each 
    x-grid cell (meters)
- `compaction_displacement_max`: Maximum compaction displacement field at each 
    x-grid cell along the pre-transport sediment-sticky interface (meters)
"""
function advect_markers_using_compaction(
    model::ModelData,
    sediment_thickness_markers_initial::Vector{Float64},
    compaction_displacement_max::Vector{Float64}
)::Nothing
    sticky_thickness_markers_initial = calculate_sticky_thickness_from_markers(model)
    
    make_debug_plots = false
    if make_debug_plots
        plot_advection_input(
            model,
            compaction_displacement_max,
            sediment_thickness_markers_initial,
            sticky_thickness_markers_initial
        )
    end
    
    advect_sediment(
        model,
        sediment_thickness_markers_initial,
        compaction_displacement_max
    )
    
    advect_sticky(
        model,
        sticky_thickness_markers_initial,
        compaction_displacement_max
    )
    
    return nothing
end

""" Advect sediment markers using compaction displacement field.

# Updated Arrays
- `marker_y`: Array of marker y-coordinates (meters)

# Arguments
- `model`: Model data
- `sediment_thickness_markers_initial`: Initial sediment thickness markers at each 
    x-grid cell (meters)
- `compaction_displacement_max`: Maximum compaction displacement field at each 
    x-grid cell along the pre-transport sediment-sticky interface (meters)
"""
function advect_sediment(
    model::ModelData,
    sediment_thickness_markers_initial::Vector{Float64},
    compaction_displacement_max::Vector{Float64}
)::Nothing
    displacement_factors = calculate_sediment_compaction_displacement_factors(
        model, sediment_thickness_markers_initial
    )
    marker_displacement = calculate_sediment_marker_displacement(
        model, displacement_factors, compaction_displacement_max
    )
    move_sediment_markers_using_compaction_field(model, marker_displacement)
    return nothing
end

""" Advect sticky markers using compaction displacement field.

# Updated Arrays
- `marker_y`: Array of marker y-coordinates (meters)

# Arguments
- `model`: Model data
- `sticky_thickness_markers_initial`: Initial sticky thickness markers at each 
    x-grid cell (meters)
- `compaction_displacement_max`: Maximum compaction displacement field at each 
    x-grid cell along the pre-transport sediment-sticky interface (meters)
"""
function advect_sticky(
    model::ModelData,
    sticky_thickness_markers_initial::Vector{Float64},
    compaction_displacement_max::Vector{Float64}
)::Nothing
    displacement_factors = calculate_sticky_compaction_displacement_factors(
        model, sticky_thickness_markers_initial
    )
    marker_displacement = calculate_sticky_marker_displacement(
        model, displacement_factors, compaction_displacement_max
    )
    move_sticky_markers_using_compaction_field(model, marker_displacement)
    return nothing
end

""" Calculate compaction displacement factors.

Displacement factors are unit vectors that indicate the displacement of
sediment markers relative to the sediment water interface due to
compaction. Vertical displacement is maximum at the
sediment-sticky-air/water interface.

# Arguments
- `model`: Model data
- `sediment_thickness_markers`: Initial sediment thickness markers at each 
    x-grid cell (meters)

# Returns
- `marker_displacement_factors`: Displacement factors for sediment markers
"""
function calculate_sediment_compaction_displacement_factors(
    model::ModelData,
    sediment_thickness_markers::Vector{Float64}
)::Vector{Float64}
    marker_arrays = model.markers.arrays
    marker_matids = marker_arrays.material.marker_matid.array
    location = marker_arrays.location
    marker_x = location.marker_x.array
    marker_y = location.marker_y.array
    gridt = model.topography.arrays.gridt.array
    gridx = gridt[1, :]
    
    matids_sedimentary_basin = get_matids_sedimentary_basin(model)
    
    use_power_law_displacement = true # Set to false for linear displacement
    
    marknum = length(marker_x)
    marker_displacement_factors = zeros(marknum)
    for imarker in 1:marknum
        matid = marker_matids[imarker]
        if matid in matids_sedimentary_basin
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            y_mudline = get_depth(x_marker, gridt)
            y_submud = y_marker - y_mudline
            sediment_thickness = linear_interp_at_x_location(
                x_marker, gridx, sediment_thickness_markers
            )
            if y_submud > sediment_thickness
                sediment_thickness = y_submud
            end
            if abs(sediment_thickness) != 0.0
                if use_power_law_displacement
                    displacement_factor = power_law_displacement(
                        y_submud, sediment_thickness
                    )
                else
                    displacement_factor = linear_displacement(
                        y_submud, sediment_thickness
                    )
                end
            else
                displacement_factor = 0.0
            end
            check_advection_model_parameters(
                y_submud, sediment_thickness, displacement_factor
            )
            marker_displacement_factors[imarker] = displacement_factor
        end
    end
    return marker_displacement_factors
end

""" Check compaction advection model parameters.

# Arguments
- `y_submud`: Distance below mudline (meters)
- `sediment_thickness`: Sediment thickness (meters)
- `displacement_factor`: Displacement factor
"""
function check_advection_model_parameters(
    y_submud::Float64,
    sediment_thickness::Float64,
    displacement_factor::Float64
)::Nothing
    if 0.0 > displacement_factor > 1.0
        println("y_submud: ", y_submud)
        println("sediment_thickness: ", sediment_thickness)
        println("displacement_factor: ", displacement_factor)
        error("Displacement factor must be between 0 and 1.")
    end
    if sediment_thickness < 0.0
        println("y_submud: ", y_submud)
        println("sediment_thickness: ", sediment_thickness)
        println("displacement_factor: ", displacement_factor)
        error("Sediment thickness must be positive.")
    end
    return nothing
end

""" Calculate linear displacement factor.

# Arguments
- `y_submud`: Distance below mudline (meters)
- `sediment_thickness`: Sediment thickness (meters)

# Returns
- `displacement_factor`: Linear displacement factor
"""
function linear_displacement(
    y_submud::Float64,
    sediment_thickness::Float64
)::Float64
    displacement_factor = 1.0 - y_submud / sediment_thickness
    return displacement_factor
end

""" Calculate power law displacement factor that approximates compaction.

The exponent scales with total sediment thickness to ensure that deep
compacted sediments are experiences limited displacement.

# Arguments
- `y_submud`: Distance below mudline (meters)
- `sediment_thickness`: Sediment thickness (meters)

# Returns
- `displacement_factor`: Power law displacement factor
"""
function power_law_displacement(
    y_submud::Float64,
    sediment_thickness::Float64
)::Float64
    exponent = max(1.0, sediment_thickness/2000.0)
    displacement_factor = (1.0 - y_submud / sediment_thickness)^exponent
    return displacement_factor
end

""" Calculate marker displacement.

# Arguments
- `model`: Model data
- `marker_displacement_factors`: Displacement factors for sediment markers
- `compaction_displacement_max`: Maximum compaction displacement field at each 
    x-grid cell along the pre-transport sediment-sticky interface (meters)

# Returns
- `marker_displacement`: Displacement of sediment markers due to compaction (meters)
"""
function calculate_sediment_marker_displacement(
    model::ModelData,
    marker_displacement_factors::Vector{Float64},
    compaction_displacement_max::Vector{Float64}
)::Vector{Float64}
    marker_arrays = model.markers.arrays
    marker_matids = marker_arrays.material.marker_matid.array
    location = marker_arrays.location
    marker_x = location.marker_x.array
    gridt = model.topography.arrays.gridt.array
    gridx = gridt[1, :]
    
    matids_sedimentary_basin = get_matids_sedimentary_basin(model)
    
    marknum = length(marker_x)
    marker_displacement = zeros(marknum)
    for imarker in 1:marknum
        matid = marker_matids[imarker]
        if matid in matids_sedimentary_basin
            x_marker = marker_x[imarker]
            displacement_factor = marker_displacement_factors[imarker]
            displacement_max = linear_interp_at_x_location(
                x_marker, gridx, compaction_displacement_max
            )
            displacement = displacement_factor * displacement_max
            marker_displacement[imarker] = displacement
        end
    end
    return marker_displacement
end

""" Calculate marker displacement.

# Updated Arrays
- `marker_y`: Array of marker y-coordinates (meters)

# Arguments
- `model`: Model data
- `marker_displacement`: Displacement of sediment markers due to compaction (meters)
"""
function move_sediment_markers_using_compaction_field(
    model::ModelData,
    marker_displacement::Vector{Float64}
)::Nothing
    marker_arrays = model.markers.arrays
    marker_matids = marker_arrays.material.marker_matid.array
    marker_y = marker_arrays.location.marker_y.array
    
    matids_sedimentary_basin = get_matids_sedimentary_basin(model)
    
    marknum = length(marker_y)
    for imarker in 1:marknum
        if marker_matids[imarker] in matids_sedimentary_basin
            marker_y[imarker] = marker_y[imarker] + marker_displacement[imarker]
        end
    end
    return nothing
end

""" Calculate sticky thickness from markers.

# Arguments
- `model`: Model data

# Returns
- `sticky_thickness`: Sticky thickness (meters)
"""
function calculate_sticky_thickness_from_markers(
    model::ModelData
)::Vector{Float64}
    gridt = model.topography.arrays.gridt.array
    gridx = gridt[1, :]
    
    matids_sticky = get_sticky_matids(model)
    material_ids_of_layer = [matids_sticky[1], matids_sticky[2]]
    
    sticky_thickness = calculate_layer_thickness(
        model, material_ids_of_layer, gridx, use_smoothing=true
    )
    return sticky_thickness
end

""" Calculate compaction displacement factors for sticky air/water.

Displacement factors are unit vectors that indicate the displacement of
sticky markers relative to the sediment water interface due to
compaction. Vertical displacement is maximum at the interface between
sticky air/water and sediment.

# Arguments
- `model`: Model data
- `sticky_thickness_markers`: Initial sticky thickness markers at each 
    x-grid cell (meters)

# Returns
- `marker_displacement_factors`: Displacement factors for sticky markers
"""
function calculate_sticky_compaction_displacement_factors(
    model::ModelData,
    sticky_thickness_markers::Vector{Float64}
)::Vector{Float64}
    marker_arrays = model.markers.arrays
    marker_matids = marker_arrays.material.marker_matid.array
    location = marker_arrays.location
    marker_x = location.marker_x.array
    marker_y = location.marker_y.array
    gridt = model.topography.arrays.gridt.array
    gridx = gridt[1, :]
    
    matids_sticky = get_sticky_matids(model)
    
    marknum = length(marker_x)
    marker_displacement_factors = zeros(marknum)
    for imarker in 1:marknum
        if marker_matids[imarker] in matids_sticky
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            sticky_thickness = linear_interp_at_x_location(
                x_marker, gridx, sticky_thickness_markers
            )
            if abs(sticky_thickness) != 0.0
                displacement_factor = y_marker / sticky_thickness
            else
                displacement_factor = 0.0
            end
            marker_displacement_factors[imarker] = displacement_factor
        end
    end
    return marker_displacement_factors
end

""" Calculate marker displacement.

# Arguments
- `model`: Model data
- `marker_displacement_factors`: Displacement factors for sticky markers
- `compaction_displacement_max`: Maximum compaction displacement field at each 
    x-grid cell along the pre-transport sediment-sticky interface (meters)

# Returns
- `marker_displacement`: Displacement of sticky markers due to compaction (meters)
"""
function calculate_sticky_marker_displacement(
    model::ModelData,
    marker_displacement_factors::Vector{Float64},
    compaction_displacement_max::Vector{Float64}
)::Vector{Float64}
    marker_arrays = model.markers.arrays
    marker_matids = marker_arrays.material.marker_matid.array
    location = marker_arrays.location
    marker_x = location.marker_x.array
    gridt = model.topography.arrays.gridt.array
    gridx = gridt[1, :]
    
    matids_sticky = get_sticky_matids(model)
    
    marknum = length(marker_x)
    marker_displacement = zeros(marknum)
    for imarker in 1:marknum
        if marker_matids[imarker] in matids_sticky
            x_marker = marker_x[imarker]
            displacement_factor = marker_displacement_factors[imarker]
            displacement_max = linear_interp_at_x_location(
                x_marker, gridx, compaction_displacement_max
            )
            displacement = displacement_factor * displacement_max
            marker_displacement[imarker] = displacement
        end
    end
    return marker_displacement
end

function move_sticky_markers_using_compaction_field(
    model::ModelData,
    marker_displacement::Vector{Float64}
)::Nothing
    marker_arrays = model.markers.arrays
    marker_matids = marker_arrays.material.marker_matid.array
    marker_y = marker_arrays.location.marker_y.array
    
    matids_sticky = get_sticky_matids(model)
    
    marknum = length(marker_y)
    for imarker in 1:marknum
        if marker_matids[imarker] in matids_sticky
            @inbounds marker_y[imarker] = marker_y[imarker] + marker_displacement[imarker]
        end
    end
    return nothing
end

function get_sticky_matids(
    model::ModelData
)::Tuple{Int16, Int16}
    matid_air = model.materials.dicts.matid_types["StickyAir"][1]
    matid_water = model.materials.dicts.matid_types["StickyWater"][1]
    matids_sticky = (matid_air, matid_water)
    return matids_sticky
end

function get_matids_sedimentary_basin(
    model::ModelData
)::Vector{Int16}
    matids_sedimentary_basin = fill(-1, 2)
    matid_types = model.materials.dicts.matid_types
    if !isempty(matid_types["Sediment"])
        matids_sedimentary_basin[1] = matid_types["Sediment"][1]
    end
    if !isempty(matid_types["SolidifiedBasalt"])
        matids_sedimentary_basin[2] = matid_types["SolidifiedBasalt"][1]
    end
    return matids_sedimentary_basin
end

""" Plot advection input.

# Arguments
- `model`: Model data
- `compaction_displacement_max`: Maximum compaction displacement field at each 
    x-grid cell along the pre-transport sediment-sticky interface (meters)
- `sediment_thickness_markers_initial`: Initial sediment thickness markers at each 
    x-grid cell (meters)
- `sticky_thickness_markers_initial`: Initial sticky thickness markers at each 
    x-grid cell (meters)
"""
function plot_advection_input(
    model::ModelData,
    compaction_displacement_max::Vector{Float64},
    sediment_thickness_markers_initial::Vector{Float64},
    sticky_thickness_markers_initial::Vector{Float64}
)::Nothing
    
    gridt = model.topography.arrays.gridt.array
    gridx = gridt[1, :]
    
    p = Plots.plot(
        gridx, compaction_displacement_max .* 100,
        label="Compaction Displacement Max*100"
    )
    Plots.plot!(
        p,
        gridx, sediment_thickness_markers_initial,
        label="Sediment Thickness Markers"
    )
    Plots.plot!(
        p,
        gridx, sticky_thickness_markers_initial,
        label="Sticky Thickness Markers"
    )
    Plots.xlabel!("x (m)")
    Plots.ylabel!("Thickness or Displacement max (m)")
    
    ntimestep = model.timestep.parameters.main_time_loop.ntimestep.value
    Plots.savefig(p, "advection_input_ntimestep_$(ntimestep).png")
    
    return nothing
end

end # module 