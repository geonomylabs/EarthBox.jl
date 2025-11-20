module PorosityUpdate

include("SedimentPorosity.jl")

import EarthBox.ModelDataContainer: ModelData
import .SedimentPorosity: calculate_rock_props

""" Calculate rock properties for all sediment markers using porosity.

# Updated Arrays
## model.markers.arrays.material
- marker_rho.array: Array{Float64,1}

## model.markers.arrays.thermal
- marker_rhocp.array: Array{Float64,1}
    - Marker density in kg/m^3.

- marker_kt.array: Array{Float64,1}
    - Marker thermal conductivity in W/m/K.
"""
function update_rock_properties_for_porosity!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    gridt = model.topography.arrays.gridt.array

    (
        conductivity_water,
        density_water,
        heat_capacity_water
    ) = get_compaction_model_parameters(model)

    matids_sediments = model.materials.dicts.matid_types["Sediment"]
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_rho = model.markers.arrays.material.marker_rho.array
    marker_rhocp = model.markers.arrays.thermal.marker_rhocp.array
    marker_kt = model.markers.arrays.thermal.marker_kt.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    mat_cp = model.materials.arrays.mat_cp.array
    marker_porosity_initial = model.markers.arrays.material.marker_porosity_initial.array
    marker_decay_depth = model.markers.arrays.material.marker_decay_depth.array
    marker_max_burial_depth = model.markers.arrays.material.marker_max_burial_depth.array
    marknum = model.markers.parameters.distribution.marknum.value

    # Use Julia's parallel processing
    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            material_id = marker_matid[imarker]
            if use_sediment_porosity(matids_sediments, material_id)
                @inbounds begin
                    x_marker = marker_x[imarker]
                    y_marker = marker_y[imarker]
                    heat_capacity = mat_cp[material_id]
                    density = marker_rho[imarker]
                    conductivity = marker_kt[imarker]
                    porosity_at_mudline = Float64(marker_porosity_initial[imarker])
                    porosity_decay_depth = Float64(marker_decay_depth[imarker])
                    max_burial = Float64(marker_max_burial_depth[imarker])
                end
                (
                    conductivity, density, rhocp
                ) = calculate_rock_props(
                    heat_capacity, conductivity, density,
                    porosity_at_mudline, porosity_decay_depth,
                    conductivity_water, density_water, heat_capacity_water,
                    x_marker, y_marker, gridt, max_burial
                )
                @inbounds begin
                    marker_rho[imarker] = density
                    marker_rhocp[imarker] = rhocp
                    marker_kt[imarker] = conductivity
                end
            end
        end
    end
    return nothing
end

@inline function get_compaction_model_parameters(
    model::ModelData
)::Tuple{Float64,Float64,Float64}
    conductivity_water = model.materials.parameters.compaction.conductivity_water.value
    density_water = model.materials.parameters.compaction.density_water.value
    heat_capacity_water = model.materials.parameters.compaction.heat_capacity_water.value
    return conductivity_water, density_water, heat_capacity_water
end

@inline function use_sediment_porosity(
    matids_sediments::Array{Int16,1},
    material_id::Int16
)::Bool
    return material_id in matids_sediments
end

end # module 