module FourLinearSegments

using EarthBox.ModelDataContainer: ModelData
import EarthBox.Arrays: ArrayUtils
using ....Layers: calc_depths_for_layers
using ..TempIniStructs
using ..PlumeUpdate: update_temperature_for_plume!

function initialize!(model::ModelData)::Nothing
    layer_thickness = LayerThickness(
        model.geometry.parameters.sticky_air_geometry.thick_air.value,
        model.geometry.parameters.earth_layering.thick_upper_crust.value,
        model.geometry.parameters.earth_layering.thick_lower_crust.value,
        model.geometry.parameters.earth_layering.thick_upper_lith.value,
        model.geometry.parameters.earth_layering.thick_middle_lith.value,
        model.geometry.parameters.earth_layering.thick_lower_lith.value,
        -9999.0  # Not used
    )

    steady_state = model.heat_equation.parameters.steady_state
    perturbation = Perturbation(
        steady_state.amplitude_perturbation.value,
        steady_state.width_perturbation.value
    )

    temperature_bcs = TemperatureBCs(
        model.bcs.parameters.temperature.temperature_top.value,
        model.bcs.parameters.temperature.temperature_bottom.value,
        steady_state.temperature_base_lith.value
    )

    internal_temperature = InternalTemperature(
        steady_state.temperature_surface.value,
        steady_state.temperature_moho.value
    )

    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_temperature = model.markers.arrays.thermal.marker_TK.array

    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value

    mxnum = model.markers.parameters.distribution.mxnum.value
    mynum = model.markers.parameters.distribution.mynum.value
    imarker = 1  # Julia uses 1-based indexing
    for _ixm in 1:mxnum
        for _iym in 1:mynum
            x_marker = marker_x[imarker]
            y_marker = marker_y[imarker]
            temperature_kelvins = calculate_temperature(
                x_marker, y_marker, xsize, ysize,
                temperature_bcs, internal_temperature,
                layer_thickness, perturbation)
            marker_temperature[imarker] = temperature_kelvins
            imarker += 1
        end
    end

    update_temperature_for_plume!(model)
    return nothing
end

""" Calculate marker temperature for a layered model with linear segments.
Use simple linear geotherm with different gradients in crust, mantle and asthenosphere.
"""
function calculate_temperature(
    x_marker::Float64,
    y_marker::Float64,
    xsize::Float64,
    ysize::Float64,
    temperature_bcs::TemperatureBCs,
    internal_temperature::InternalTemperature,
    layer_thickness::LayerThickness,
    perturbation::Perturbation
)::Float64
    temperature_top = temperature_bcs.temperature_top
    temperature_surface = internal_temperature.temperature_surface
    temperature_moho = internal_temperature.temperature_moho
    temperature_base_lith = temperature_bcs.temperature_base_lith
    temperature_bottom = temperature_bcs.temperature_bottom

    (
        depth_layer1, _depth_layer2, depth_layer3,
        _depth_layer4, _depth_layer5, depth_layer6
    ) = calc_depths_for_layers(
        layer_thickness.thick_air,
        layer_thickness.thick_upper_crust,
        layer_thickness.thick_lower_crust,
        layer_thickness.thick_upper_lith,
        layer_thickness.thick_middle_lith,
        layer_thickness.thick_lower_lith
    )

    depth_perturbation = calculate_lithosphere_depth_perturbation(
        perturbation.amplitude_perturbation,
        perturbation.width_perturbation,
        xsize,
        x_marker
    )

    depth_layer6 = depth_layer6 - depth_perturbation

    temperature_kelvins = temperature_top
    # Air
    if y_marker <= depth_layer1
        temperature_kelvins = temperature_top
    # Crust
    elseif depth_layer3 >= y_marker > depth_layer1
        temperature_kelvins = (
            temperature_surface
            + (temperature_moho - temperature_surface)
            / (depth_layer3 - depth_layer1)
            * (y_marker - depth_layer1)
        )
    # Mantle lithosphere
    elseif depth_layer6 >= y_marker > depth_layer3
        temperature_kelvins = (
            temperature_moho
            + (temperature_base_lith - temperature_moho)
            / (depth_layer6 - depth_layer3)
            * (y_marker - depth_layer3)
        )
    # Asthenosphere
    elseif y_marker > depth_layer6
        temperature_kelvins = (
            temperature_base_lith
            + (temperature_bottom - temperature_base_lith)
            / (ysize - depth_layer6)
            * (y_marker - depth_layer6)
        )
    end
    return temperature_kelvins
end

function calculate_lithosphere_depth_perturbation(
    amplitude::Float64,
    width::Float64,
    xsize::Float64,
    x_marker::Float64
)::Float64
    xmid = xsize/2.0
    xstart = xmid - width
    xfinal = xmid + width
    depth_perturbation = 0.0
    if xfinal >= x_marker >= xstart
        if x_marker <= xmid
            depth_perturbation = amplitude/(xmid - xstart)*(x_marker - xstart)
        else
            depth_perturbation = (
                amplitude - amplitude/(xfinal - xmid)*(x_marker - xmid)
            )
        end
    end
    return depth_perturbation
end

end # module 