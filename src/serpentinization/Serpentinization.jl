module Serpentinization

include("core/UpdateRockProps.jl")
include("core/SerpentineRheology.jl")

import SpecialFunctions: erf
import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit, print_info
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.Markers.MarkerMaterials: MaterialGroupIDs
import EarthBox.SurfaceProcesses.MarkerTransformation: transform_marker_to_serpentinite!
import EarthBox.ConversionFuncs: celsius_to_kelvin
import EarthBox.SurfaceProcesses: calculate_age_ma
import EarthBox: GridFuncs
import EarthBox.SurfaceProcesses: TopoInterpolation
import .SerpentineRheology: update_marker_flow_viscosity!
import .SerpentineRheology: update_marker_failure_properties!
import .UpdateRockProps: update_rock_props!

const PDATA = get_eb_parameters()

struct ValidInputNames
    iuse_serpentinization::Symbol
    serpentinization_temperature::Symbol
    maximum_serpentinization_depth::Symbol
    maximum_serpentinization_rate::Symbol
    nominal_strain_rate_serpentinization::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize parameters for the serpentinization model.

# Arguments
- `model::ModelData`: Model data container containing the model parameters and arrays.

# Keyword Arguments
- `$(PDATA.iuse_serpentinization.name)::Int64`
    - $(PDATA.iuse_serpentinization.description)

- `$(PDATA.serpentinization_temperature.name)::Float64`
    - $(PDATA.serpentinization_temperature.description)

- `$(PDATA.maximum_serpentinization_depth.name)::Float64`
    - $(PDATA.maximum_serpentinization_depth.description)

- `$(PDATA.maximum_serpentinization_rate.name)::Float64`
    - $(PDATA.maximum_serpentinization_rate.description)

- `$(PDATA.nominal_strain_rate_serpentinization.name)::Float64`
    - $(PDATA.nominal_strain_rate_serpentinization.description)

"""
function initialize!(model::ModelData; kwargs...)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

""" Apply serpentinization model.

This method manages functions that perform the following tasks:

1. Calculations of marker serpentinization ratio.
2. Update of exothermic heat production from serpentinization.
"""
function update_marker_serpentinization!(model::ModelData, inside_flags::Vector{Int8})
    iuse_serpentinization = model.materials.parameters.serpentinization.iuse_serpentinization.value
    if iuse_serpentinization == 1
        @timeit_memit "Finished updating marker serpentinization" begin
            marker_incremental_serpentinization_ratio = calculate_marker_serpentinization(
                model, inside_flags)
            update_exothermic_heat_production_from_serpentinization!(
                model, marker_incremental_serpentinization_ratio)
        end
    end
end

function update_rock_props_for_serpentinization!(model::ModelData, inside_flags::Vector{Int8})
        marker_serpentinization = model.markers.arrays.material.marker_serpentinization.array
        max_serpentinization = maximum(marker_serpentinization)
        mantle_ids = get_mantle_serpentinization_ids(model)
        iuse_serpentinization = model.materials.parameters.serpentinization.iuse_serpentinization.value
        if iuse_serpentinization == 1
            if max_serpentinization > 0.0
                @timeit_memit "Finished updating marker rock properties for serpentinization" begin
                    update_rock_props!(model, mantle_ids, inside_flags)
                end
            else
                print_info("No serpentinization detected", level=2)
            end
        end
end

function update_marker_failure_properties!(model::ModelData)::Nothing
    marker_serpentinization = model.markers.arrays.material.marker_serpentinization.array
    max_serpentinization = maximum(marker_serpentinization)
    mantle_ids = get_mantle_serpentinization_ids(model)
    if _use_serpentinization_rheology(model)
        if max_serpentinization > 0.0
            update_marker_failure_properties!(model, mantle_ids)
        end
    end
    return nothing
end

function update_marker_flow_viscosity!(model::ModelData)::Nothing
    marker_serpentinization = model.markers.arrays.material.marker_serpentinization.array
    max_serpentinization = maximum(marker_serpentinization)
    mantle_ids = get_mantle_serpentinization_ids(model)
    if _use_serpentinization_rheology(model)
        if max_serpentinization > 0.0
            update_marker_flow_viscosity!(model, mantle_ids)
        end
    end
end

function _use_serpentinization_rheology(model::ModelData)::Bool
    iuse_serpentinization = model.materials.parameters.serpentinization.iuse_serpentinization.value
    iuse_serpentinization_rheology = 0
    return iuse_serpentinization == 1 && iuse_serpentinization_rheology == 1
end

""" Calculate marker serpentinization.

See Merzi et al. (2024) for details.

# Arguments
- `model`: ModelData instance

# Returns
- `marker_incremental_serpentinization_ratio`: Incremental serpentinization ratio 
  for each marker
"""
function calculate_marker_serpentinization(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Vector{Float64}
    timestep = model.timestep.parameters.main_time_loop.timestep.value

    serpentinization = model.materials.parameters.serpentinization
    serpentinization_temperature_celsius = serpentinization.serpentinization_temperature.value
    maximum_serpentinization_depth = serpentinization.maximum_serpentinization_depth.value
    maximum_serpentinization_rate = serpentinization.maximum_serpentinization_rate.value
    nominal_strain_rate_serpentinization = serpentinization.nominal_strain_rate_serpentinization.value

    serpentinization_temperature_kelvins = celsius_to_kelvin(serpentinization_temperature_celsius)

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_TK = model.markers.arrays.thermal.marker_TK.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marker_serpentinization = model.markers.arrays.material.marker_serpentinization.array
    marker_strain_rate_plastic = model.markers.arrays.strain.marker_strain_rate_plastic.array

    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    dx_topo = model.topography.parameters.topo_grid.dx_topo.value
    toponum = model.topography.parameters.topo_grid.toponum.value
    gridt = model.topography.arrays.gridt.array

    mantle_ids = get_mantle_serpentinization_ids(model)
    marknum = model.markers.parameters.distribution.marknum.value

    marker_incremental_serpentinization_ratio = Vector{Float64}(undef, marknum)

    Threads.@threads for imarker in 1:marknum
        incremental_serpentinization_ratio = 0.0
        if inside_flags[imarker] == 1
            @inbounds begin
                matid = marker_matid[imarker]
                y_marker = marker_y[imarker]
            end
            if matid in mantle_ids && y_marker >= y_sealevel
                @inbounds x_marker = marker_x[imarker]
                submud_depth = calculate_submud_depth(
                    x_marker, y_marker, gridt, dx_topo, toponum)
                if submud_depth < maximum_serpentinization_depth
                    @inbounds temperature_kelvins = marker_TK[imarker]
                    strain_rate_plastic = marker_strain_rate_plastic[imarker]
                    serpentinization_ratio_initial = marker_serpentinization[imarker]
                    incremental_serpentinization_ratio = 
                        calculate_incremental_serpentinization_ratio(
                            serpentinization_ratio_initial, 
                            strain_rate_plastic,
                            temperature_kelvins, 
                            timestep,
                            serpentinization_temperature_kelvins,
                            maximum_serpentinization_rate,
                            nominal_strain_rate_serpentinization
                        )
                    
                    marker_serpentinization[imarker] = serpentinization_ratio_initial + incremental_serpentinization_ratio

                end
            end
        end
        @inbounds marker_incremental_serpentinization_ratio[imarker] = 
            incremental_serpentinization_ratio
    end

    min_serp = minimum(marker_serpentinization)
    max_serp = maximum(marker_serpentinization)
    min_incr_serp = minimum(marker_incremental_serpentinization_ratio)
    max_incr_serp = maximum(marker_incremental_serpentinization_ratio)

    print_info("Min marker_serpentinization: $(min_serp)", level=2)
    print_info("Max marker_serpentinization: $(max_serp)", level=2)
    print_info("Min marker_incremental_serpentinization_ratio: $(min_incr_serp)", level=2)
    print_info("Max marker_incremental_serpentinization_ratio: $(max_incr_serp)", level=2)

    return marker_incremental_serpentinization_ratio
end

""" Calculate incremental serpentinization ratio.

See Merzi et al. (2024) for details.

# Arguments
- `serpentinization_ratio_initial`: Initial serpentinization ratio
- `strain_rate_plastic`: Plastic strain rate (1/s)
- `temperature_kelvins`: Temperature in Kelvin
- `timestep`: Timestep in seconds
- `serpentinization_temperature_kelvins`: Serpentinization temperature in Kelvin 
  where reaction rate reaches maximum and then rapidly decreases
- `maximum_serpentinization_rate`: Maximum serpentinization rate
- `nominal_strain_rate_serpentinization`: Nominal plastic strain rate for 
  serpentinization

# Returns
- `incremental_serpentinization_ratio`: Incremental serpentinization ratio
"""
@inline function calculate_incremental_serpentinization_ratio(
    serpentinization_ratio_initial::Float64,
    strain_rate_plastic::Float64,
    temperature_kelvins::Float64,
    timestep::Float64,
    serpentinization_temperature_kelvins::Float64,
    maximum_serpentinization_rate::Float64,
    nominal_strain_rate_serpentinization::Float64
)::Float64
    
    serpentinization_rate = calculate_serpentinization_rate(
        temperature_kelvins,
        strain_rate_plastic,
        maximum_serpentinization_rate,
        serpentinization_temperature_kelvins,
        nominal_strain_rate_serpentinization
    )

    ratio_new = serpentinization_rate * (1.0 - serpentinization_ratio_initial) * timestep
    ratio_limit = 1.0 - serpentinization_ratio_initial

    return min(ratio_new, ratio_limit)
end

""" Calculate serpentinization rate.

See Merzi et al. (2024) for details.
"""
@inline function calculate_serpentinization_rate(
    temperature_kelvins::Float64,
    strain_rate_plastic::Float64,
    maximum_serpentinization_rate::Float64,
    serpentinization_temperature_kelvins::Float64,
    nominal_strain_rate_serpentinization::Float64
)::Float64
    
    temperature_dependent_rate_factor = calculate_temperature_dependent_rate_factor(
        temperature_kelvins, serpentinization_temperature_kelvins)

    plastic_strain_rate_rate_factor = calculate_plastic_strain_rate_rate_factor(
        strain_rate_plastic, nominal_strain_rate_serpentinization)

    return maximum_serpentinization_rate * 
           temperature_dependent_rate_factor * 
           plastic_strain_rate_rate_factor
end

""" Calculate temperature dependent rate factor.

See Merzi et al. (2024) for details.

# Arguments
- `temperature_kelvins`: Temperature in Kelvin
- `temperature_reference_kelvins`: Reference temperature in Kelvin

# Returns
- `temperature_dependent_rate_factor`: Temperature dependent rate factor
"""
@inline function calculate_temperature_dependent_rate_factor(
    temperature_kelvins::Float64,
    temperature_reference_kelvins::Float64
)::Float64
    
    alpha_factor = -4.0
    w_factor = 94.0
    x_factor = (temperature_kelvins - temperature_reference_kelvins) / w_factor
    theta_factor = exp(-(x_factor^2.0) / 2.0)
    phi_factor = 0.5 * (1.0 + erf(alpha_factor * x_factor / sqrt(2.0)))
    
    return 2.0 * theta_factor * phi_factor
end

""" Calculate plastic strain rate rate factor.

See Merzi et al. (2024) for details.

# Arguments
- `strain_rate_plastic`: Plastic strain rate (1/s)
- `nominal_strain_rate_serpentinization`: Reference strain rate (1/s)

# Returns
- `plastic_strain_rate_rate_factor`: Plastic strain rate rate factor
"""
@inline function calculate_plastic_strain_rate_rate_factor(
    strain_rate_plastic::Float64,
    nominal_strain_rate_serpentinization::Float64
)::Float64
    strain_rate_ratio = strain_rate_plastic / nominal_strain_rate_serpentinization
    return 1.0 - exp(-strain_rate_ratio)
end

""" Get all material IDs of mantle melting materials.

# Arguments
- `model`: ModelData instance
- `use_serpentinite`: Whether to include serpentinite IDs

# Returns
- `mantle_melting_mat_ids`: Material IDs of mantle melting materials
"""
function get_mantle_serpentinization_ids(
    model::ModelData;
    use_serpentinite::Bool=false
)::Vector{Int16}
    matid_types = model.materials.dicts.matid_types
    fertile = matid_types["UltramaficMantleFertile"]
    nfertile = length(fertile)
    partially_molten = matid_types["UltramaficMantlePartiallyMolten"]
    npartially_molten = length(partially_molten)
    refractory = matid_types["UltramaficMantleRefactory"]
    nrefractory = length(refractory)
    ntotal = nfertile + npartially_molten + nrefractory
    if use_serpentinite
        serpentinite = matid_types["Serpentinite"]
        nserpentinite = length(serpentinite)
        ntotal += nserpentinite
    end
    mantle_melting_mat_ids = zeros(Int16, ntotal)
    icount = 1
    for i in 1:nfertile
        mantle_melting_mat_ids[icount] = fertile[i]
        icount += 1
    end
    for i in 1:npartially_molten
        mantle_melting_mat_ids[icount] = partially_molten[i]
        icount += 1
    end
    for i in 1:nrefractory
        mantle_melting_mat_ids[icount] = refractory[i]
        icount += 1
    end
    if use_serpentinite
        for i in 1:nserpentinite
            mantle_melting_mat_ids[icount] = serpentinite[i]
            icount += 1
        end
    end
    return mantle_melting_mat_ids
end

""" Calculate exothermic heat production from serpentinization.

The enthalpy change for the exothermic transformation of olivine (forsterite)
into serpentine (chrysotile) is -19.5 kcal/mol (e.g. Merzi et al., 2024).
"""
function update_exothermic_heat_production_from_serpentinization!(
    model::ModelData,
    marker_incremental_serpentinization_ratio::Vector{Float64}
)::Nothing
    enthalpy_kcal_mol = -19.5 # kcal/mol
    enthalpy_joules_mol = enthalpy_kcal_mol * 4184.0 # J/mol
    molar_volume_m3_mol = get_molar_volume_of_chrysotile()

    mxstep = model.markers.parameters.distribution.mxstep.value
    mystep = model.markers.parameters.distribution.mystep.value
    avg_marker_volume = mxstep * mystep

    marknum = length(marker_incremental_serpentinization_ratio)
    timestep = model.timestep.parameters.main_time_loop.timestep.value

    marker_serpentinization_heat_production = 
        model.markers.arrays.material.marker_serpentinization_heat_production.array

    Threads.@threads for imarker in 1:marknum
        @inbounds incremental_serpentinization_ratio = marker_incremental_serpentinization_ratio[imarker]
        incremental_serpentine_volume = avg_marker_volume * incremental_serpentinization_ratio
        incremental_moles = incremental_serpentine_volume / molar_volume_m3_mol
        # Use abs so negative exothermic leads to positive heat production
        total_energy = abs(incremental_moles * enthalpy_joules_mol)
        power = total_energy / timestep # Watts
        power_density = power / avg_marker_volume # W/m^3
        @inbounds marker_serpentinization_heat_production[imarker] = power_density
    end

    print_serpentinization_info(
        marker_incremental_serpentinization_ratio,
        marker_serpentinization_heat_production
    )
    return nothing
end

function print_serpentinization_info(
    marker_incremental_serpentinization_ratio::Vector{Float64},
    marker_serpentinization_heat_production::Vector{Float64}
)::Nothing
    min_serpentinization_ratio = minimum(marker_incremental_serpentinization_ratio)
    max_serpentinization_ratio = maximum(marker_incremental_serpentinization_ratio)
    print_info("Min incremental serpentinization ratio: $(min_serpentinization_ratio)", level=2)
    print_info("Max incremental serpentinization ratio: $(max_serpentinization_ratio)", level=2)
    
    min_heat_production = minimum(marker_serpentinization_heat_production)
    max_heat_production = maximum(marker_serpentinization_heat_production)
    print_info("Min serpentinization heat production (W/m^3): $(min_heat_production)", level=2)
    print_info("Max serpentinization heat production (W/m^3): $(max_heat_production)", level=2)
    return nothing
end

""" Calculate molar volume of chrysotile.

# Returns
- `molar_volume_m3_mol`: Molar volume of chrysotile in m^3/mol
"""
function get_molar_volume_of_chrysotile()::Float64
    molar_mass_g_mol = get_molar_mass_of_chrysotile()
    molar_mass_kg_mol = molar_mass_g_mol / 1000.0
    density_kg_m3 = 2500.0 # kg/m^3
    return molar_mass_kg_mol / density_kg_m3
end

""" Calculate molar mass of chrysotile.

The chemical formula for chrysotile is Mg3Si2O5(OH)4:

- Magnesium (Mg): 3 atoms x 24.305 g/mol = 72.915 g/mol
- Silicon (Si): 2 atoms x 28.085 g/mol = 56.170 g/mol
- Oxygen (O): 9 atoms x 16.00 g/mol = 144.00 g/mol
- Hydrogen (H): 4 atoms x 1.008 g/mol = 4.032 g/mol

molar_mass = 72.915 + 56.170 + 144.00 + 4.032 = 277.117 g/mol

# Returns
- `molar_mass_chrysotile`: Molar mass of chrysotile in g/mol
"""
@inline function get_molar_mass_of_chrysotile()::Float64
    return 277.117 # g/mol
end

@inline function calculate_submud_depth(
    x_marker::Float64,
    y_marker::Float64,
    gridt::Matrix{Float64},
    dx_topo::Float64,
    toponum::Int
)::Float64
    ixn, dx_dist = TopoInterpolation.find_closest_topography_node_to_marker(
        x_marker, gridt, dx_topo, toponum)
    y_topo = TopoInterpolation.calculate_topography(dx_dist, gridt, ixn)
    return y_marker - y_topo
end

end # module 