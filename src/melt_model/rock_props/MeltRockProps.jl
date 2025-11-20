module MeltRockProps

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs
import ..MeltPropertiesOpt: marker_melt_props

""" Calculate rock properties for all markers.

# Updated Arrays
## model.markers.arrays.material
- marker_rho.array::Vector{Float64}
    - Marker density in kg/m^3.

## model.markers.arrays.thermal
- marker_rhocp.array::Vector{Float64}
    - Marker density in kg/m^3.

- marker_ha.array::Vector{Float64}
    - Adiabatic term (expansivity x temperature) for marker.

# Returns
- marker_props::Matrix{Float64}: (4, marknum)
    - Marker properties:
        - temperature_kelvins
        - pressure_pascals
        - delta_heat_capacity
        - delta_expansivity

"""
function update_for_melt!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    iuse_melt_thermal_props = model.melting.parameters.options.iuse_melt_thermal_props.value
    iuse_depletion_density = model.melting.parameters.options.iuse_depletion_density.value
    marknum = model.markers.parameters.distribution.marknum.value
    # Marker quantity arrays
    marker_TK = model.markers.arrays.thermal.marker_TK.array
    marker_pr = model.markers.arrays.pressure.marker_pr.array
    marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array
    marker_extracted_meltfrac = model.markers.arrays.melt.marker_extracted_meltfrac.array
    marker_extractable_meltfrac = model.markers.arrays.melt.marker_extractable_meltfrac.array
    marker_rho = model.markers.arrays.material.marker_rho.array
    marker_rhocp = model.markers.arrays.thermal.marker_rhocp.array
    marker_ha = model.markers.arrays.thermal.marker_ha.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    # Material Arrays
    mat_melting_itypes = model.materials.arrays.mat_melting_itypes.array
    mat_melting = model.materials.arrays.mat_melting.array
    mat_rho = model.materials.arrays.mat_rho.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                # Marker quantities
                pressure_pascals = marker_pr[imarker]
                temperature_kelvins = marker_TK[imarker]
                _melt_fraction = marker_meltfrac[imarker]
                extracted_melt_fraction = marker_extracted_meltfrac[imarker]
                extractable_melt_fraction = marker_extractable_meltfrac[imarker]
                density = marker_rho[imarker]
                rhocp = marker_rhocp[imarker]
                adiabatic_term = marker_ha[imarker]

                # Material parameters associated with marker
                matid = marker_matid[imarker]
                density_melt = mat_rho[matid, 4]
                itype_solidus = mat_melting_itypes[matid, 1]
                itype_liquidus = mat_melting_itypes[matid, 2]
                latent_heat = mat_melting[matid, 1]
            end

            (
                density, rhocp, adiabatic_term, 
                delta_heat_capacity, delta_expansivity
            ) = marker_melt_props(
                    pressure_pascals,
                    temperature_kelvins,
                    density,
                    rhocp,
                    adiabatic_term,
                    iuse_melt_thermal_props,
                    density_melt,
                    _melt_fraction,
                    extracted_melt_fraction,
                    extractable_melt_fraction,
                    itype_solidus,
                    itype_liquidus,
                    latent_heat,
                    iuse_depletion_density
                )
            @inbounds begin
                marker_rho[imarker] = density
                marker_rhocp[imarker] = rhocp
                marker_ha[imarker] = adiabatic_term
            end
        end
    end
    return nothing
end

end # module 