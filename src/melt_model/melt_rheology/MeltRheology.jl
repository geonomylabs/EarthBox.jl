module MeltRheology

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Markers.MarkerMaterials.MaterialGroupIDs: get_sticky_material_ids
import EarthBox: GridFuncs
import ..Extraction: get_mantle_melting_ids

function update_marker_melt_rheology!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    iuse_exponential_viscosity_reduction = 
        model.melting.parameters.options.iuse_exponential_viscosity_reduction.value
    # First calculate the effect of partial melt on the viscosity of the rocks
    if iuse_exponential_viscosity_reduction == 1
        println("   >> Using exponential viscosity reduction for partially molten rocks.")
        update_marker_melt_rheology_exponential(model, inside_flags)
    else
        println("   >> Setting viscosity to melt viscosity if melt fraction exceeds 10%.")
        update_marker_melt_rheology_using_threshold(model, inside_flags)
    end

    # Next calculate the effects of melt damage due to transport from partially
    # molten domain to the emplacement and extrusive domains
    iuse_melt_damage = model.materials.parameters.melt_damage.iuse_melt_damage.value
    if iuse_melt_damage == 1
        println("   >> Applying melt damage to marker flow viscosity.")
        update_marker_flow_viscosity_for_melt_damage(model, inside_flags)
    end
    return nothing
end

""" Update marker viscosity and strain based on melt fraction using a threshold approach.

# Arguments
- `model::ModelData`: Model data container object

# Updated Arrays
- `model.markers.arrays.rheology.marker_eta`: Marker viscosity (Pa.s)
- `model.markers.arrays.strain.marker_GI`: Accumulated strain in markers
- `model.markers.arrays.strain.marker_strain_plastic`: Accumulated plastic strain
- `model.markers.arrays.strain.marker_sr_ratio`: Ratio of strain rate calculated 
  using grid stress changes and a Maxwell model over strain rate interpolated 
  from the grid
"""
function update_marker_melt_rheology_using_threshold(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    marknum = model.markers.parameters.distribution.marknum.value
    marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array
    marker_sr_ratio = model.markers.arrays.strain.marker_sr_ratio.array
    marker_GII = model.markers.arrays.strain.marker_GII.array
    marker_strain_plastic = model.markers.arrays.strain.marker_strain_plastic.array
    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_extracted_meltfrac = model.markers.arrays.melt.marker_extracted_meltfrac.array
    viscosity_melt = model.melting.parameters.rheology.viscosity_melt.value

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                mextracted_meltfrac = marker_extracted_meltfrac[imarker]
                xmelt = marker_meltfrac[imarker]
            end
            # Only use melting model if incremental melt extraction >= 0.
            # Otherwise the material is considered refractory.
            if xmelt > 0 && timesum > 0 && mextracted_meltfrac >= 0
                # Viscosity of partially molten rocks
                if xmelt > 0.1
                    @inbounds marker_eta[imarker] = viscosity_melt
                end
                @inbounds begin
                    # Reset strain rate ratio
                    marker_sr_ratio[imarker] = 1
                    # Reset strain
                    marker_GII[imarker] = 0
                    marker_strain_plastic[imarker] = 0
                end
            end
        end
    end
    return nothing
end

""" Update marker viscosity and strain based on melt fraction using an exponential approach.

# Arguments
- `model::ModelData`: Model data container object
- `alpha_factor::Float64`: Factor used in exponential term for reducing effective 
  viscosity of partially molten rocks
- `xmelt_threshold::Float64`: Threshold melt fraction above which effective 
  viscosity is set to viscosity_melt
"""
function update_marker_melt_rheology_exponential(
    model::ModelData,
    inside_flags::Vector{Int8};
    alpha_factor::Float64=35.0,
    xmelt_threshold::Float64=0.3
)
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    marknum = model.markers.parameters.distribution.marknum.value
    viscosity_melt = model.melting.parameters.rheology.viscosity_melt.value

    mantle_melting_ids = get_mantle_melting_ids(model)

    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array
    marker_extractable_meltfrac = model.markers.arrays.melt.marker_extractable_meltfrac.array
    marker_extracted_meltfrac = model.markers.arrays.melt.marker_extracted_meltfrac.array
    marker_sr_ratio = model.markers.arrays.strain.marker_sr_ratio.array
    marker_GII = model.markers.arrays.strain.marker_GII.array
    marker_strain_plastic = model.markers.arrays.strain.marker_strain_plastic.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                matid = marker_matid[imarker]
                xmelt = marker_meltfrac[imarker]
                xmelt_extracted = marker_extracted_meltfrac[imarker]
                xmelt_extractable = marker_extractable_meltfrac[imarker]
            end
            # Only use melting model if incremental melt extraction >= 0.
            # Otherwise the material is considered refractory.
            if xmelt > 0 && timesum > 0 && xmelt_extracted >= 0
                # Mantle rocks with melt extraction
                if matid in mantle_melting_ids
                    if xmelt_extractable > 0.0
                        # Use xmelt_extractable since extraction occurs in
                        # mantle rocks and only the extractable melt fraction
                        # is present in the rock
                        @inbounds begin
                            marker_eta[imarker] = update_viscosity_with_partial_melt(
                                marker_eta[imarker], xmelt_extractable,
                                alpha_factor, xmelt_threshold, viscosity_melt
                                )
                            # Reset strain rate ratio
                            marker_sr_ratio[imarker] = 1
                            # Reset strain
                            marker_GII[imarker] = 0
                            marker_strain_plastic[imarker] = 0
                        end
                    end
                # Crustal rocks without melt extraction
                else
                    # Use xmelt since no extraction occurs in crustal rocks
                    @inbounds begin
                        marker_eta[imarker] = update_viscosity_with_partial_melt(
                            marker_eta[imarker], xmelt, alpha_factor,
                            xmelt_threshold, viscosity_melt
                            )
                        # Reset strain rate ratio
                        marker_sr_ratio[imarker] = 1
                        # Reset strain
                        marker_GII[imarker] = 0
                        marker_strain_plastic[imarker] = 0
                    end
                end
            end
        end
    end
end

""" Update the effective viscosity of partially molten rocks.

# Arguments
- `effective_viscosity::Float64`: Marker viscosity (Pa.s)
- `xmelt::Float64`: Melt fraction
- `alpha_factor::Float64`: Factor used in exponential term for reducing effective 
  viscosity
- `xmelt_threshold::Float64`: Threshold melt fraction above which effective 
  viscosity is set to viscosity_melt
- `viscosity_melt::Float64`: Viscosity of melt

# Returns
- `Float64`: Updated effective viscosity
"""
@inline function update_viscosity_with_partial_melt(
    effective_viscosity::Float64,
    xmelt::Float64,
    alpha_factor::Float64,
    xmelt_threshold::Float64,
    viscosity_melt::Float64
)::Float64
    if xmelt < xmelt_threshold
        effective_viscosity = exponential_melt_rheology(
            effective_viscosity, xmelt, alpha_factor)
    else
        effective_viscosity = viscosity_melt
    end
    return effective_viscosity
end

""" Calculate the effective viscosity of partially molten rocks using an exponential term.

# Arguments
- `marker_eta::Float64`: Marker viscosity (Pa.s)
- `xmelt::Float64`: Melt fraction
- `alpha_factor::Float64`: Factor used in exponential term for reducing effective 
  viscosity

# Returns
- `Float64`: Effective viscosity of partially molten rocks
"""
@inline function exponential_melt_rheology(
    marker_eta::Float64,
    xmelt::Float64,
    alpha_factor::Float64
)::Float64
    return marker_eta / exp(alpha_factor * xmelt)
end

""" Update marker viscosity and strain based on melt fraction for melt damage.

# Arguments
- `model::ModelData`: Model data container object
- `reset_strain::Bool`: Whether to reset strain values

# Updated Arrays
- `model.markers.arrays.rheology.marker_eta`: Marker viscosity (Pa.s)
- `model.markers.arrays.strain.marker_GI`: Accumulated strain in markers
- `model.markers.arrays.strain.marker_strain_plastic`: Accumulated plastic strain
- `model.markers.arrays.strain.marker_sr_ratio`: Ratio of strain rate calculated 
  using grid stress changes and a Maxwell model over strain rate interpolated 
  from the grid
"""
function update_marker_flow_viscosity_for_melt_damage(
    model::ModelData,
    inside_flags::Vector{Int8};
    reset_strain::Bool=false
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    viscosity_melt = model.melting.parameters.rheology.viscosity_melt.value

    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_sr_ratio = model.markers.arrays.strain.marker_sr_ratio.array
    marker_GII = model.markers.arrays.strain.marker_GII.array
    marker_strain_plastic = model.markers.arrays.strain.marker_strain_plastic.array
    marker_melt_damage = model.markers.arrays.strain.marker_melt_damage.array
    sticky_ids = get_sticky_material_ids(model)

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                matid = marker_matid[imarker]
                marker_damage_factor = marker_melt_damage[imarker]
            end
            if matid ∉ sticky_ids && marker_damage_factor > 1
                @inbounds begin
                    eta_damaged = marker_eta[imarker] / marker_damage_factor
                    marker_eta[imarker] = max(eta_damaged, viscosity_melt)
                end
                if reset_strain
                    @inbounds begin
                        # Reset strain rate ratio
                        marker_sr_ratio[imarker] = 1
                        # Reset strain
                        marker_GII[imarker] = 0
                        marker_strain_plastic[imarker] = 0
                    end
                end
            end
        end
    end
    return nothing
end

end # module 