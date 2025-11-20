module SedimentThickness

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelStructureManager.TopAndBottom: calculate_layer_thickness

""" Calculate sediment thickness from markers including extrusives.

# Arguments
- `model`: Model data
- `topo_gridx`: Topography grid x-coordinates (meters)

# Returns
- `sediment_thickness`: Sediment thickness (meters)
"""
function calculate_sediment_thickness_from_markers(
    model::ModelData,
    topo_gridx::Vector{Float64}
)::Vector{Float64}
    matid_types = model.materials.dicts.matid_types
    matid_sediment = matid_types["Sediment"][1]
    
    matid_basalt = Int16(-1)
    if !isempty(matid_types["SolidifiedBasalt"])
        matid_basalt = matid_types["SolidifiedBasalt"][1]
    end
    
    matid_salt = Int16(-1)
    if !isempty(matid_types["Salt"])
        matid_salt = matid_types["Salt"][1]
    end
    
    material_ids_of_layer = [matid_sediment, matid_basalt, matid_salt]
    sediment_thickness = calculate_layer_thickness(
        model, material_ids_of_layer, topo_gridx, use_smoothing=true
    )
    return sediment_thickness
end

""" Calculate sediment thickness from markers including extrusives.

# Arguments
- `model`: Model data

# Returns
- `sediment_thickness`: Sediment thickness (meters)
"""
function calculate_sediment_thickness_from_markers_no_basalt(
    model::ModelData
)::Vector{Float64}
    gridt = model.topography.arrays.gridt.array
    gridx = gridt[1, :]
    matid_sediment = model.materials.dicts.matid_types["Sediment"][1]
    material_ids_of_layer = [matid_sediment]
    sediment_thickness = calculate_layer_thickness(
        model, material_ids_of_layer, gridx, use_smoothing=true
    )
    return sediment_thickness
end

end # module