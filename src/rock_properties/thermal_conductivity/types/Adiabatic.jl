module Adiabatic

using EarthBox.ModelDataContainer: ModelData

"""
    update_conductivity!(model::ModelData)

Update thermal conductivity for all markers using the adiabatic formulation.

# Arguments
- `model::ModelData`: The model data container containing marker information

# Updates
- `model.markers.arrays.thermal.marker_kt.array`: Array of marker thermal conductivities
"""
function update_conductivity!(model::ModelData, inside_flags::Vector{Int8})
    marker_kt = model.markers.arrays.thermal.marker_kt.array
    marker_matid = model.markers.arrays.material.marker_matid.array
    marknum = model.markers.parameters.distribution.marknum.value

    Threads.@threads for imarker in 1:marknum
        if inside_flags[imarker] == 1
            @inbounds begin
                material_id = marker_matid[imarker]
                y_location = model.markers.arrays.location.marker_y.array[imarker]
                conductivity_reference = model.materials.arrays.mat_kt.array[material_id, 1]
            end
            conductivity = marker_conductivity_model(conductivity_reference, y_location)
            @inbounds marker_kt[imarker] = conductivity
        end
    end
end

"""
    marker_conductivity_model(
        conductivity_reference::Float64,
        y_location::Float64
    )::Float64

Calculate marker thermal conductivity using the adiabatic formulation.

# Arguments
- `conductivity_reference::Float64`: Reference thermal conductivity
- `y_location::Float64`: Depth of marker in meters

# Returns
- `Float64`: Calculated thermal conductivity
"""
@inline function marker_conductivity_model(
    conductivity_reference::Float64,
    y_location::Float64
)::Float64
    return thermal_conductivity_simple_lithosphere_100km(y_location, conductivity_reference)
end

"""
    thermal_conductivity_simple_lithosphere_100km(
        y_location::Float64,
        thermal_conductivity::Float64
    )::Float64

Force adiabatic gradient below lithosphere using thermal conductivity.

# Arguments
- `y_location::Float64`: Depth of marker below top of model in meters
- `thermal_conductivity::Float64`: Thermal conductivity (W/m/K)

# Returns
- `Float64`: Thermal conductivity (W/m/K)
"""
@inline function thermal_conductivity_simple_lithosphere_100km(
    y_location::Float64,
    thermal_conductivity::Float64
)::Float64
    lithosphere_thickness = 100000.0
    sticky_air_thickness = 10000.0
    lithosphere_asthenosphere_boundary_depth = lithosphere_thickness + sticky_air_thickness

    if y_location > lithosphere_asthenosphere_boundary_depth
        asthenosphere_thickness = 10000.0
        temperature_top = 0.0
        temperature_base = 1334.5
        temperature_base_lithosphere = 1330.0
        lithosphere_thermal_conductivity = 2.25

        thermal_conductivity = calc_thermal_conductivity_for_forced_adiabate(
            lithosphere_thickness,
            asthenosphere_thickness,
            temperature_top,
            temperature_base,
            temperature_base_lithosphere,
            lithosphere_thermal_conductivity
        )
    end

    return thermal_conductivity
end

"""
    calc_thermal_conductivity_for_forced_adiabate(
        lithosphere_thickness::Float64,
        asthenosphere_thickness::Float64,
        temperature_top::Float64,
        temperature_base::Float64,
        temperature_base_lithosphere::Float64,
        lithosphere_thermal_conductivity::Float64
    )::Float64

Calculate thermal conductivity required to maintain adiabatic gradient.

# Arguments
- `lithosphere_thickness::Float64`: Thickness of lithosphere in meters
- `asthenosphere_thickness::Float64`: Thickness of asthenosphere in meters
- `temperature_top::Float64`: Temperature at the top of the lithosphere in Celsius
- `temperature_base::Float64`: Temperature at the base of the model domain in Celsius
- `temperature_base_lithosphere::Float64`: Temperature at the base of the lithosphere in Celsius
- `lithosphere_thermal_conductivity::Float64`: Thermal conductivity of the lithosphere in W/m/K

# Returns
- `Float64`: Thermal conductivity in W/m/K in the sub-lithospheric domain required
             to produce the adiabatic gradient assuming only conductive heat transport

# Background
Assuming a lithosphere with thickness L, a sub-lithospheric domain with
thickness A, a lithosphere thermal conductivity of kL, temperature at
the top of the lithosphere Ttop, temperature at the base of the model
Tbase and temperature at the base of the lithopshere TBL, the
thermal conductivity kA in the sublithospheric domain required to
produce a constant steady-state geotherm in the sub-lithospheric domain
is given by:

    kA = A/L*kL*(TBL - Ttop)/(Tb - TBL)
"""
@inline function calc_thermal_conductivity_for_forced_adiabate(
    lithosphere_thickness::Float64,
    asthenosphere_thickness::Float64,
    temperature_top::Float64,
    temperature_base::Float64,
    temperature_base_lithosphere::Float64,
    lithosphere_thermal_conductivity::Float64
)::Float64
    return (asthenosphere_thickness / lithosphere_thickness) *
           lithosphere_thermal_conductivity *
           (temperature_base_lithosphere - temperature_top) /
           (temperature_base - temperature_base_lithosphere)
end

end # module 