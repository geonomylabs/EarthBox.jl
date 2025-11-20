module MarkersStrain

import EarthBox.ModelDataContainer: ModelData
import EarthBox: GridFuncs

""" Update marker total strain.

# Updated arrays from group `markers.arrays.strain`
- `marker_GII::Vector{Float64}`: Marker strain
"""
function update_total_strain!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    marknum = model.markers.parameters.distribution.marknum.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    marker_exx = model.markers.arrays.strain.marker_exx.array
    marker_exy = model.markers.arrays.strain.marker_exy.array
    marker_strain = model.markers.arrays.strain.marker_GII.array
    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                exx = marker_exx[imarker]
                exy = marker_exy[imarker]
                strain = marker_strain[imarker]
                marker_strain[imarker] = calculate_strain(exx, exy, timestep, strain)
            end
        end
    end
    return nothing
end

"""
    update_plastic_strain!(model::ModelData)

Update marker plastic strain.

# Updated arrays from group `markers.arrays.strain`
- `marker_strain_plastic::Vector{Float64}`: Marker strain
- `marker_strain_rate_plastic::Vector{Float64}`: Marker strain rate
"""
function update_plastic_strain!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    marker_sxx = model.markers.arrays.stress.marker_sxx.array
    marker_sxy = model.markers.arrays.stress.marker_sxy.array
    marker_eta = model.markers.arrays.rheology.marker_eta.array
    marker_eta_flow = model.markers.arrays.rheology.marker_eta_flow.array
    marker_strain_plastic = model.markers.arrays.strain.marker_strain_plastic.array
    marker_strain_rate_plastic = model.markers.arrays.strain.marker_strain_rate_plastic.array
    marker_pfailure = model.markers.arrays.rheology.marker_pfailure.array

    (
        plastic_healing_rate
    ) = model.materials.parameters.stress_limits_yield.plastic_healing_rate.value

    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            if marker_pfailure[imarker] == 1
                viscoplastic_viscosity = marker_eta[imarker]
                flow_viscosity = marker_eta_flow[imarker]
                sxx = marker_sxx[imarker]
                sxy = marker_sxy[imarker]
                plastic_strain_old = marker_strain_plastic[imarker]
                
                strain_plastic, strain_rate_plastic = calculate_plastic_strain(
                    plastic_strain_old,
                    timestep,
                    sxx,
                    sxy,
                    viscoplastic_viscosity,
                    flow_viscosity,
                    plastic_healing_rate
                )
                
                marker_strain_plastic[imarker] = strain_plastic
                marker_strain_rate_plastic[imarker] = strain_rate_plastic
            end
        end
    end
    return nothing
end

@inline function calculate_strain(
    strain_rate_xx::Float64,
    strain_rate_xy::Float64,
    timestep::Float64,
    strain_old::Float64
)::Float64
    strain_new = strain_old + timestep * sqrt(strain_rate_xx^2 + strain_rate_xy^2)
    return strain_new
end

""" Calculate plastic strain based on stress and viscosities.

# Updated arrays from group `markers.arrays.strain`
- `marker_strain_plastic::Vector{Float64}`: Marker strain
- `marker_strain_rate_plastic::Vector{Float64}`: Marker strain rate

# Arguments
- `plastic_strain_old::Float64`: Plastic strain at previous timestep
- `timestep::Float64`: Timestep (seconds)
- `sxx::Float64`: Normal stress (Pa)
- `sxy::Float64`: Shear stress (Pa)
- `viscoplastic_viscosity::Float64`: Viscoplastic viscosity (Pa.s)
- `flow_viscosity::Float64`: Flow viscosity (Pa.s)
- `plastic_healing_rate::Float64`: Plastic healing rate (1/s)

# Returns
- `Tuple{Float64, Float64}`: (plastic_strain, plastic_strain_rate)
"""
function calculate_plastic_strain(
    plastic_strain_old::Float64,
    timestep::Float64,
    sxx::Float64,
    sxy::Float64,
    viscoplastic_viscosity::Float64,
    flow_viscosity::Float64,
    plastic_healing_rate::Float64
)::Tuple{Float64, Float64}
    stress_invariant = calculate_stress_invariant(sxx, sxy)
    plastic_strain_rate = calculate_plastic_strain_rate(
        stress_invariant,
        viscoplastic_viscosity,
        flow_viscosity
    )
    plastic_strain_rate = max(0.0, plastic_strain_rate)
    plastic_strain = plastic_strain_old + timestep * plastic_strain_rate - 
                    plastic_healing_rate * timestep
    plastic_strain = max(0.0, plastic_strain)
    return plastic_strain, plastic_strain_rate
end

function calculate_stress_invariant(sxx::Float64, sxy::Float64)::Float64
    stress_invariant = sqrt(sxx^2 + sxy^2)
    return stress_invariant
end

function calculate_plastic_strain_rate(
    stress_invariant::Float64,
    viscoplastic_viscosity::Float64,
    flow_viscosity::Float64
)::Float64
    plastic_strain_rate = 0.5 * stress_invariant * (
        1.0 / viscoplastic_viscosity - 1.0 / flow_viscosity
    )
    return plastic_strain_rate
end

end # module 