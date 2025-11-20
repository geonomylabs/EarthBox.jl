module AnalyticalThreeLayer

using EarthBox.ModelDataContainer: ModelData
using ....Layers: calc_depths_for_layers
using ..Geotherm
using ..TempIniStructs
using ..PlumeUpdate: update_temperature_for_plume!
using ..GeothermStructs: SteadyState3LayerData

function initialize!(model::ModelData)::Nothing

    geometry_params = model.geometry.parameters
    thick_air = geometry_params.sticky_air_geometry.thick_air.value
    thick_upper_crust = geometry_params.earth_layering.thick_upper_crust.value
    thick_lower_crust = geometry_params.earth_layering.thick_lower_crust.value

    @assert !isnan(thick_air) "Thickness of sticky air (`thick_air`) is not set"
    @assert !isnan(thick_upper_crust) "Thickness of upper crust (`thick_upper_crust`) is not set"
    @assert !isnan(thick_lower_crust) "Thickness of lower crust (`thick_lower_crust`) is not set"

    steady_state_params = model.heat_equation.parameters.steady_state
    
    # Do not use asserts since defaults are set to useful values
    conductivity_upper_crust = steady_state_params.conductivity_upper_crust.value
    conductivity_lower_crust = steady_state_params.conductivity_lower_crust.value
    conductivity_mantle = steady_state_params.conductivity_mantle.value
    heat_production_upper_crust = steady_state_params.heat_production_upper_crust.value
    heat_production_lower_crust = steady_state_params.heat_production_lower_crust.value
    heat_production_mantle = steady_state_params.heat_production_mantle.value
    amplitude_perturbation = steady_state_params.amplitude_perturbation.value
    width_perturbation = steady_state_params.width_perturbation.value

    thick_thermal_lithosphere = steady_state_params.thick_thermal_lithosphere.value
    temperature_base_lith = steady_state_params.temperature_base_lith.value
    @assert !isnan(thick_thermal_lithosphere) "Thickness of thermal lithosphere (`thick_thermal_lithosphere`) is not set"
    @assert !isnan(temperature_base_lith) "Temperature at base of lithosphere (`temperature_base_lith`) is not set"

    bc_temperature_params = model.bcs.parameters.temperature
    temperature_top = bc_temperature_params.temperature_top.value
    temperature_bottom = bc_temperature_params.temperature_bottom.value # Not used by the AnalyticalThreeLayer model

    adiabatic_gradient = model.heat_equation.parameters.initial_condition.adiabatic_gradient.value

    layer_thickness = LayerThickness(
        thick_air,
        thick_upper_crust,
        thick_lower_crust,
        0.0, # thick_upper_lith is not used by the AnalyticalThreeLayer model
        0.0, # thick_middle_lith is not used by the AnalyticalThreeLayer model
        0.0, # thick_lower_lith is not used by the AnalyticalThreeLayer model
        thick_thermal_lithosphere
    )

    conductivity = Conductivity(
        conductivity_upper_crust,
        conductivity_lower_crust,
        conductivity_mantle
    )

    heat_production = HeatProduction(
        heat_production_upper_crust,
        heat_production_lower_crust,
        heat_production_mantle
    )

    temperature_bcs = TemperatureBCs(
        temperature_top,
        0.0, # temperature_bottom is not used by the AnalyticalThreeLayer model
        temperature_base_lith
    )

    perturbation = Perturbation(
        amplitude_perturbation,
        width_perturbation
    )

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_temperature = model.markers.arrays.thermal.marker_TK.array

    xsize = model.grids.parameters.geometry.xsize.value

    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    imarker = 1
    for _ixm in 1:mxnum
        for _iym in 1:mynum
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            temperature_kelvins = calculate_temperature(
                x_marker, y_marker, xsize, layer_thickness,
                conductivity, heat_production, temperature_bcs,
                perturbation, adiabatic_gradient)
            marker_temperature[imarker] = temperature_kelvins
            imarker += 1
        end
    end

    update_temperature_for_plume!(model)
    return nothing
end

"""
    calculate_temperature(
        x_marker::Float64, 
        y_marker::Float64, 
        xsize::Float64,
        layer_thickness::LayerThickness, 
        conductivity::Conductivity,
        heat_production::HeatProduction, 
        temperature_bcs::TemperatureBCs,
        perturbation::Perturbation, 
        adiabatic_gradient::Float64
    )::Float64

Apply perturbations and calculate marker temperature.
"""
function calculate_temperature(
    x_marker::Float64,
    y_marker::Float64,
    xsize::Float64,
    layer_thickness::LayerThickness,
    conductivity::Conductivity,
    heat_production::HeatProduction,
    temperature_bcs::TemperatureBCs,
    perturbation::Perturbation,
    adiabatic_gradient::Float64
)::Float64
    (
        depth_air,
        depth_asthenosphere_adjusted,
        thick_thermal_lithosphere_adjusted,
        _thick_thermal_mantle_lithosphere_adjusted
    ) = calc_geometry_for_3layer_model(
        x_marker, xsize, layer_thickness, perturbation)

    temperature_kelvins = calc_steadystate_temperature_for_3layer_model(
        y_marker,
        depth_air,
        depth_asthenosphere_adjusted,
        thick_thermal_lithosphere_adjusted,
        layer_thickness,
        conductivity,
        heat_production,
        temperature_bcs,
        adiabatic_gradient
    )
    return temperature_kelvins
end

"""
    calc_geometry_for_3layer_model(
        x_marker::Float64, xsize::Float64,
        layer_thickness::LayerThickness,
        perturbation::Perturbation
    )::Tuple{Float64,Float64,Float64,Float64}

Calculate the geometry of the 3-layer lithosphere model taking into account 
perturbations to the thermal lithosphere.
"""
@inline function calc_geometry_for_3layer_model(
    x_marker::Float64,
    xsize::Float64,
    layer_thickness::LayerThickness,
    perturbation::Perturbation
)::Tuple{Float64,Float64,Float64,Float64}
    (
        depth_air, _depth_layer2, _depth_layer3,
        _depth_layer4, _depth_layer5, _depth_layer6
    ) = calc_depths_for_layers(
        layer_thickness.thick_air,
        layer_thickness.thick_upper_crust,
        layer_thickness.thick_lower_crust,
        layer_thickness.thick_upper_lith,
        layer_thickness.thick_middle_lith,
        layer_thickness.thick_lower_lith
    )

    depth_asthenosphere = layer_thickness.thick_thermal_lithosphere + depth_air

    thick_thermal_mantle_lithosphere = (
        layer_thickness.thick_thermal_lithosphere -
        layer_thickness.thick_upper_crust -
        layer_thickness.thick_lower_crust
    )

    (
        depth_asthenosphere_adjusted,
        thick_thermal_lithosphere_adjusted,
        thick_thermal_mantle_lithosphere_adjusted
    ) = calculate_adjusted_thermal_lithosphere_geometry(
        x_marker,
        depth_asthenosphere,
        thick_thermal_mantle_lithosphere,
        xsize,
        layer_thickness.thick_thermal_lithosphere,
        perturbation
    )

    return (
        depth_air,
        depth_asthenosphere_adjusted,
        thick_thermal_lithosphere_adjusted,
        thick_thermal_mantle_lithosphere_adjusted
    )
end

"""
    calculate_adjusted_thermal_lithosphere_geometry(
        x_marker::Float64,
        depth_asthenosphere::Float64,
        thick_thermal_mantle_lithosphere::Float64,
        xsize::Float64,
        thick_thermal_lithosphere::Float64,
        perturbation::Perturbation
    )::Tuple{Float64,Float64,Float64}

Calculate adjusted lithosphere geometry.
"""
function calculate_adjusted_thermal_lithosphere_geometry(
    x_marker::Float64,
    depth_asthenosphere::Float64,
    thick_thermal_mantle_lithosphere::Float64,
    xsize::Float64,
    thick_thermal_lithosphere::Float64,
    perturbation::Perturbation
)::Tuple{Float64,Float64,Float64}
    amplitude = perturbation.amplitude_perturbation
    width = perturbation.width_perturbation
    xmid = xsize/2.0
    xss = xmid - width
    xff = xmid + width
    telev = 0.0
    if xss <= x_marker <= xff
        if x_marker <= xmid
            telev = amplitude/(xmid - xss)*(x_marker - xss)
        else
            telev = amplitude - amplitude/(xff - xmid)*(x_marker - xmid)
        end
    end
    depth_asthenosphere_adjusted = depth_asthenosphere - telev
    thick_thermal_lithosphere_adjusted = thick_thermal_lithosphere - telev
    thick_thermal_mantle_lithosphere_adjusted = thick_thermal_mantle_lithosphere - telev
    return (
        depth_asthenosphere_adjusted,
        thick_thermal_lithosphere_adjusted,
        thick_thermal_mantle_lithosphere_adjusted
    )
end

"""
    calc_steadystate_temperature_for_3layer_model(
        y_marker::Float64,
        depth_air::Float64,
        depth_asthenosphere::Float64,
        thick_thermal_lithosphere::Float64,
        layer_thickness::LayerThickness,
        conductivity::Conductivity,
        heat_production::HeatProduction,
        temperature_bcs::TemperatureBCs,
        adiabatic_grad::Float64;
        use_adiabatic_gradient::Bool=true
    )::Float64

Calculate steady-state temperature for 3-layer lithosphere.
"""
@inline function calc_steadystate_temperature_for_3layer_model(
    y_marker::Float64,
    depth_air::Float64,
    depth_asthenosphere::Float64,
    thick_thermal_lithosphere::Float64,
    layer_thickness::LayerThickness,
    conductivity::Conductivity,
    heat_production::HeatProduction,
    temperature_bcs::TemperatureBCs,
    adiabatic_grad::Float64;
    use_adiabatic_gradient::Bool=true
)::Float64
    thick_thermal_mantle_lithosphere = (
        thick_thermal_lithosphere -
        layer_thickness.thick_upper_crust -
        layer_thickness.thick_lower_crust
    )

    temperature_kelvins = temperature_bcs.temperature_top
    if y_marker <= depth_air
        temperature_kelvins = temperature_bcs.temperature_top
    elseif depth_air < y_marker <= depth_asthenosphere
        z_depth_sub_rock = y_marker - depth_air

        ss_3layer_data = SteadyState3LayerData(
            temperature_bcs.temperature_top,
            temperature_bcs.temperature_base_lith,
            layer_thickness.thick_upper_crust,
            layer_thickness.thick_lower_crust,
            thick_thermal_mantle_lithosphere,
            thick_thermal_lithosphere,
            heat_production.heat_production_upper_crust,
            heat_production.heat_production_lower_crust,
            heat_production.heat_production_mantle,
            conductivity.conductivity_upper_crust,
            conductivity.conductivity_lower_crust,
            conductivity.conductivity_mantle
        )
        (
            temperature_kelvins,
            _surface_heat_flow
        ) = Geotherm.steady_state_3layer(z_depth_sub_rock, ss_3layer_data)
    else
        if use_adiabatic_gradient
            delta_temperature_adiabatic = (y_marker - depth_asthenosphere)*adiabatic_grad/1000.0
        else
            delta_temperature_adiabatic = 0.0
        end
        temperature_kelvins = temperature_bcs.temperature_base_lith + delta_temperature_adiabatic
    end
    return temperature_kelvins
end

end # module 