module GetModelData

import EarthBox.ModelDataContainer: ModelData
import ..DataNames: MarkerDataNames
import ..DataNames: get_list
import ..ArrayLookup: create_array_lookup, get_array_info
import ...PlotParametersManager: PlotParameters
import ...PlotParametersManager.PlotConversionManager: convert_time_units
import ...PlotParametersManager.PlotConversionManager: convert_length_units
import ...PlotParametersManager.PlotConversionManager: convert_length_array_units

function get_model_marker_data(
    parameters::PlotParameters,
    marker_data_names::MarkerDataNames,
    model_data::ModelData
)::Tuple{Float64, Float64, Float64, Dict{String, Vector{Float64}}}
    
    model_time = model_data.timestep.parameters.main_time_loop.timesum.value
    time_units = model_data.timestep.parameters.main_time_loop.timesum.units

    y_sealevel_meters = model_data.topography.parameters.sealevel.y_sealevel.value
    y_sealevel = convert_length_units(parameters.conversion, 'm', y_sealevel_meters)

    base_level_shift_meters = model_data.topography.parameters.sealevel.base_level_shift.value
    base_level_shift = convert_length_units(parameters.conversion, 'm', base_level_shift_meters)

    (
        model_time, time_units
    ) = convert_time_units(parameters.conversion, model_time, time_units)

    array_lookup = create_array_lookup(model_data)
    marker_data_dict = Dict{String, Vector{Float64}}()
    
    for dataname in get_list(marker_data_names)
        array, units = get_array_info(array_lookup, dataname)
        if dataname in ["marker_x", "marker_y"]
            array = convert_length_array_units(parameters.conversion, units, array)
        end
        marker_data_dict[dataname] = array
    end
    
    return model_time, y_sealevel, base_level_shift, marker_data_dict
end

end # module
