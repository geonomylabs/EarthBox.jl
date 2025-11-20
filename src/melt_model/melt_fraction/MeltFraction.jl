module MeltFraction

import EarthBox.ModelDataContainer: ModelData
import EarthBox.GridFuncs: check_in_domain
import EarthBox.SolidusLiquidus: Sediment
import EarthBox.SolidusLiquidus: ContinentalCrustUpper
import EarthBox.SolidusLiquidus: ContinentalCrustLower
import EarthBox.SolidusLiquidus: Mantle
import EarthBox.SolidusLiquidus: Gabbro

""" Update melt fraction.

Melt fraction is a function of marker temperature, pressure and material type.

# Updated Array
- model.markers.arrays.melt.marker_meltfrac.array::Vector{Float64} (marknum)
"""
function update_marker_melt_fraction!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    marker_pr = model.markers.arrays.pressure.marker_pr.array
    marker_temperature = model.markers.arrays.thermal.marker_TK.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array
    mat_melting_itypes = model.materials.arrays.mat_melting_itypes.array
    mat_melting = model.materials.arrays.mat_melting.array

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                matid = marker_matid[imarker]
                itype_solidus = mat_melting_itypes[matid, 1]
                itype_liquidus = mat_melting_itypes[matid, 2]
                latent_heat = mat_melting[matid, 1]
                pr = marker_pr[imarker]
                temperature = marker_temperature[imarker]
            end
            melt_fraction, _latent_heat_increment = 
                calc_melt_fraction_and_melt_fraction_times_latent_heat(
                    pr,
                    temperature,
                    itype_solidus,
                    itype_liquidus,
                    latent_heat
                )
            @inbounds marker_meltfrac[imarker] = melt_fraction
        end
    end
    return nothing
end

""" Calculate melt fraction.

This function computes melt fraction and melt fraction times latent heat at
given pressure and temperature.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals
- `temperature_kelvins::Float64`: Marker temperature in Kelvins
- `itype_solidus::Int`: Solidus type id
- `itype_liquidus::Int`: Liquidus type id
- `latent_heat::Float64`: Latent heat in J/kg

# Returns
- `melt_fraction::Float64`: Melt fraction of marker
- `latent_heat_increment::Float64`: Latent heat increment in J/kg
"""
@inline function calc_melt_fraction_and_melt_fraction_times_latent_heat(
    pressure_pascals::Float64,
    temperature_kelvins::Float64,
    itype_solidus::Int,
    itype_liquidus::Int,
    latent_heat::Float64
)::Tuple{Float64, Float64}
    temperature_liquidus, temperature_solidus = get_melting_model_parameters(
        pressure_pascals,
        itype_solidus,
        itype_liquidus
    )
    # default values for melt fraction and incremental latent heat
    melt_fraction = 0.0
    latent_heat_increment = 0.0
    # A liquidus equal to -1 indicates that no melting model was defined and
    # that the marker will not melt
    if temperature_liquidus > 0
        temperature_solidus = prevent_solidus_liquidus_intersection(
            temperature_solidus, temperature_liquidus)
        melt_fraction = calculate_melt_fraction(
            temperature_kelvins,
            temperature_solidus,
            temperature_liquidus
        )
        latent_heat_increment = calc_latent_heat_increment(
            latent_heat, melt_fraction)
    end
    return melt_fraction, latent_heat_increment
end

""" Prevent solidus and liquidus intersection.

# Arguments
- `temperature_solidus::Float64`: Solidus temperature in Kelvins
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins

# Returns
- `temperature_solidus::Float64`: Solidus temperature in Kelvins
"""
@inline function prevent_solidus_liquidus_intersection(
    temperature_solidus::Float64,
    temperature_liquidus::Float64
)::Float64
    # Solidus and liquidus must not intersect in the extrapolation region
    if temperature_solidus > temperature_liquidus - 100.0
        temperature_solidus = temperature_liquidus - 100.0
    end
    return temperature_solidus
end

""" Calculate melt fraction and apply limits

# Arguments
- `temperature_kelvins::Float64`: Marker temperature in Kelvins
- `temperature_solidus::Float64`: Solidus temperature in Kelvins
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins

# Returns
- `melt_fraction::Float64`: Melt fraction
"""
@inline function calculate_melt_fraction(
    temperature_kelvins::Float64,
    temperature_solidus::Float64,
    temperature_liquidus::Float64
)::Float64
    melt_fraction = (
        (temperature_kelvins - temperature_solidus)
        / (temperature_liquidus - temperature_solidus)
    )
    melt_fraction = clamp(melt_fraction, 0.0, 1.0)
    return melt_fraction
end

""" Calculate latent heat increment.

This function computes latent heat increment (latent heat x melt fraction).

# Arguments
- `latent_heat::Float64`: Latent heat in J/kg
- `melt_fraction::Float64`: Melt fraction of marker

# Returns
- `latent_heat_increment::Float64`: Latent heat increment in J/kg
"""
@inline function calc_latent_heat_increment(
    latent_heat::Float64,
    melt_fraction::Float64
)::Float64
    latent_heat_increment = latent_heat * melt_fraction
    return latent_heat_increment
end

""" Calculate liquidus and solidus temperatures and define latent heat.

If no melting model is defined, the function will return a liquidus equal
-1 indicating that the marker will not melt.

Equations used to calculate parameters are based on the material id.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals
- `itype_solidus::Int`: Solidus type id
- `itype_liquidus::Int`: Liquidus type id

# Returns
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins
- `temperature_solidus::Float64`: Solidus temperature in Kelvins
"""
@inline function get_melting_model_parameters(
    pressure_pascals::Float64,
    itype_solidus::Int,
    itype_liquidus::Int
)::Tuple{Float64, Float64}
    temperature_liquidus = get_liquidus(pressure_pascals, itype_liquidus)
    temperature_solidus = get_solidus(pressure_pascals, itype_solidus)
    return temperature_liquidus, temperature_solidus
end

""" Calculate solidus temperature.

A returned solidus temperature of -1 indicates a model has not been defined
and the marker will not melt.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals
- `itype_solidus::Int`: Solidus type id

# Returns
- `temperature_solidus::Float64`: Solidus temperature in Kelvins. 
    Negative one indicates no melting
"""
@inline function get_solidus(
    pressure_pascals::Float64,
    itype_solidus::Int
)::Float64
    # Use default parameters that will prevent melting
    temperature_solidus = -1.0
    # Sediment melting model
    if itype_solidus == 0
        temperature_solidus = Sediment.solidus_gerya2010(pressure_pascals)
    # Upper continental crust melting model
    elseif itype_solidus == 1
        temperature_solidus = ContinentalCrustUpper.solidus_gerya2010(pressure_pascals)
    # Lower continental crust melting model
    elseif itype_solidus == 2
        temperature_solidus = ContinentalCrustLower.solidus_gerya2010(pressure_pascals)
    # Ultramafic mantle models
    elseif itype_solidus == 3
        temperature_solidus = Mantle.solidus_gerya2010(pressure_pascals)
    elseif itype_solidus == 4
        temperature_solidus = Mantle.solidus_katz2003(pressure_pascals)
    # Gabbro melting models
    # Dry Gabbro
    elseif itype_solidus == 5
        temperature_solidus = Gabbro.solidus_gerya2010(pressure_pascals)
    # Wet Gabbro
    elseif itype_solidus == 6
        temperature_solidus = Gabbro.solidus_schubert2013(pressure_pascals)
    # Layered (fractionated) gabbro in lower oceanic crust
    elseif itype_solidus == 7
        temperature_solidus = Gabbro.solidus_gerya2010_gabbro_glacier(pressure_pascals)
    end
    return temperature_solidus
end

""" Calculate liquidus.

A returned liquidus temperature of -1 indicates a model has not been defined
and the marker will not melt.

# Arguments
- `pressure_pascals::Float64`: Marker pressure in Pascals
- `itype_liquidus::Int`: Liquidus type id

# Returns
- `temperature_liquidus::Float64`: Liquidus temperature in Kelvins. 
    Negative one indicates no melting
"""
@inline function get_liquidus(
    pressure_pascals::Float64,
    itype_liquidus::Int
)::Float64
    # Use default parameters that will prevent melting
    temperature_liquidus = -1.0
    # Sediment melting model
    if itype_liquidus == 0
        temperature_liquidus = Sediment.liquidus_gerya2010(pressure_pascals)
    # Upper continental crust melting model
    elseif itype_liquidus == 1
        temperature_liquidus = ContinentalCrustUpper.liquidus_gerya2010(pressure_pascals)
    # Lower continental crust melting model
    elseif itype_liquidus == 2
        temperature_liquidus = ContinentalCrustLower.liquidus_gerya2010(pressure_pascals)
    # Ultramafic mantle models
    elseif itype_liquidus == 3
        temperature_liquidus = Mantle.liquidus_gerya2010(pressure_pascals)
    elseif itype_liquidus == 4
        temperature_liquidus = Mantle.liquidus_katz2003(pressure_pascals)
    # Gabbro melting model
    elseif itype_liquidus == 5
        temperature_liquidus = Gabbro.liquidus_gerya2010(pressure_pascals)
    # Layered (fractionated) gabbro in lower oceanic crust
    elseif itype_liquidus == 6
        temperature_liquidus = Gabbro.liquidus_gerya2010_gabbro_glacier(pressure_pascals)
    end
    return temperature_liquidus
end

""" Get material id arrays.

This function returns material id arrays for sediment, upper continental
crust, lower continental crust, ultramafic mantle and gabbro.

# Arguments
- `model::ModelData`: ModelData instance

# Returns
- `melt_model_matids::Matrix{Int16}`: Material id array for melting model.
    melt_model_matids[1,:]: material id's for sediment
    melt_model_matids[2,:]: material id's for upper continental crust
    melt_model_matids[3,:]: material id's for lower continental crust
    melt_model_matids[4,:]: material id's for ultramafic mantle
    melt_model_matids[5,:]: material id's for gabbro
    where -1 denotes elements that are not used.
"""
@inline function get_matid_array_for_melting_model(model::ModelData)::Matrix{Int16}
    matid_types = model.materials.dicts.matid_types
    # Sediment melting model matid's
    sediment_matids = matid_types["Sediment"]
    ns = length(sediment_matids)
    # Upper continental crust modeling model matid's
    upper_continental_crust_matids = vcat(
        matid_types["FelsicContinentalCrustFertile"],
        matid_types["FelsicContinentalCrustPartiallyMolten"],
        matid_types["FelsicContinentalCrustRefactory"],
        matid_types["ExtractedFelsicMagma"],
        matid_types["SolidifiedGranite"],
        matid_types["SolidifiedRhyolite"]
    )
    nuc = length(upper_continental_crust_matids)
    # Lower continental crust melting model matid's
    lower_continental_crust_matids = vcat(
        matid_types["MaficContinentalCrustFertile"],
        matid_types["MaficContinentalCrustPartiallyMolten"],
        matid_types["MaficContinentalCrustRefactory"]
    )
    nlc = length(lower_continental_crust_matids)
    # Mantle melting model matid's
    ultramafic_mantle_matids = vcat(
        matid_types["UltramaficMantleFertile"],
        matid_types["UltramaficMantlePartiallyMolten"],
        matid_types["UltramaficMantleRefactory"]
    )
    num = length(ultramafic_mantle_matids)
    # Gabbro melting model matid's
    gabbro_matids = vcat(
        matid_types["ExtractedGabbroicMagma"],
        matid_types["ExtrudedGabbroicMagma"],
        matid_types["SolidifiedGabbro"],
        matid_types["SolidifiedGabbroPartiallyMolten"],
        matid_types["SolidifiedBasalt"]
    )
    ng = length(gabbro_matids)

    tmp_array = [ns, nuc, nlc, num, ng]
    nmax = maximum(tmp_array)

    sediment_matids = extend_arrays_to_max_with_negative_ones(
        nmax, sediment_matids)
    upper_continental_crust_matids = extend_arrays_to_max_with_negative_ones(
        nmax, upper_continental_crust_matids)
    lower_continental_crust_matids = extend_arrays_to_max_with_negative_ones(
        nmax, lower_continental_crust_matids)
    ultramafic_mantle_matids = extend_arrays_to_max_with_negative_ones(
        nmax, ultramafic_mantle_matids)
    gabbro_matids = extend_arrays_to_max_with_negative_ones(
        nmax, gabbro_matids)

    melt_model_matids = vcat(
        sediment_matids',
        upper_continental_crust_matids',
        lower_continental_crust_matids',
        ultramafic_mantle_matids',
        gabbro_matids'
    )

    return melt_model_matids
end

""" Extend arrays to maximum length with negative ones.

The purpose of this is to make all arrays the same length so they can be
stacked into a single array.

# Arguments
- `nmax::Int`: Maximum length of arrays
- `array::Vector{Int16}`: Array to extend

# Returns
- `array::Vector{Int16}`: Extended array
"""
@inline function extend_arrays_to_max_with_negative_ones(
    nmax::Int,
    array::Vector{Int16}
)::Vector{Int16}
    nelem = length(array)
    nextension = nmax - nelem
    if nextension > 0
        ext = fill(-1, nextension)
        array = vcat(array, ext)
    end
    return array
end

end # module 