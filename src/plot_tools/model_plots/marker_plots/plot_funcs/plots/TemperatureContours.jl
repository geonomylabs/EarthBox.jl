module TemperatureContours

import CairoMakie
import ....PlotParametersManager: PlotParameters
import ....GridPlotsManager.PlotScalarArraysManager: PlotScalarArrays
import ....PlotDtypes: AxesType
import ....PlotParametersManager.PlotContoursManager: PlotContours
import ....PlotParametersManager.PlotContoursManager: activate_contours!
import ....PlotParametersManager.PlotContoursManager: update_contour_levels!
import ....PlotParametersManager.PlotContoursManager: update_linewidths!
import ....PlotParametersManager.PlotContoursManager: plot_contours!
import ....PlotParametersManager.PlotContoursManager: add_to_contour_description!

# Start here
function plot_temperature_contours!(
    parameters::PlotParameters,
    scalar_arrays::PlotScalarArrays,
    axes::AxesType
)::Nothing
    marker_plot_params = parameters.marker_plot_params
    contour_interval = marker_plot_params.temperature_contour_interval
    plot_contours = marker_plot_params.plot_temperature_contours
    number_format = marker_plot_params.temperature_number_format
    rightside_up = marker_plot_params.temperature_label_rightside_up
    label = "Temperature (C)"
    
    if plot_contours == 1
        println(">> Plotting temperature contours")
        activate_contours!(parameters.contours)
        parameters.contours.iplot_contour_labels = marker_plot_params.plot_contour_labels
        parameters.contours.contour_interval = contour_interval
        parameters.contours.excluded_vals = [-10.0]
        value_min = marker_plot_params.temperature_min
        value_max = marker_plot_params.temperature_max
        
        update_contour_levels!(parameters.contours, value_min, value_max)
        update_linewidths!(
            parameters.contours, value_min, value_max,
            marker_plot_params.contour_line_width
        )
        
        color = Symbol(marker_plot_params.temperature_contour_color)
       
        labelsize = parameters.fonts.contour_label_fontsize
        plot_contours!(
            axes, 
            parameters.contours, 
            scalar_arrays.gridx, 
            scalar_arrays.gridy,
            scalar_arrays.scalar;
            color=color, 
            labelsize=labelsize,
            number_format=number_format
        )
        
        color_str = marker_plot_params.temperature_contour_color
        description = "    $color_str : $label : CI=$(contour_interval)"
        add_to_contour_description!(parameters.contours, description)
    end
    
    return nothing
end

end # module
