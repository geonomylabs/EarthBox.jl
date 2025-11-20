module PlotConversionManager

import EarthBox.UnitConversion: UnitConversionData
import EarthBox.UnitConversion: get_conversion_func
import ...PlotDtypes: PlotDictType, AbstractPlotParameterGroup
import ..PlotUnitsManager: PlotUnits

mutable struct PlotConversions <: AbstractPlotParameterGroup
    unit_conversion_data::UnitConversionData
    plot_units::PlotUnits
    function PlotConversions(plot_dict::PlotDictType)
        new(UnitConversionData(), PlotUnits(plot_dict))
    end
end

function convert_time_units(
    plot_conversions::PlotConversions,
    model_time::Float64,
    units_start::String
)::Tuple{Float64, String}
    (
        conversion_func_time, units_end
    ) = get_conversion_func_time(plot_conversions, units_start)
    model_time = conversion_func_time(model_time)
    return model_time, units_end
end

function get_conversion_func_time(
    plot_conversions::PlotConversions,
    units_start::String
)::Tuple{Function, String}
    units_end = plot_conversions.plot_units.time_units
    conv_func, units_start, units_end = get_conversion_func(
        units_start, units_end, plot_conversions.unit_conversion_data, "time_units")
    return conv_func, units_end
end

function convert_grid_arrays_to_plot_units(
    plot_conversions::PlotConversions,
    gridx::Vector{Float64},
    gridy::Vector{Float64},
    length_units::String
)::Tuple{Vector{Float64}, Vector{Float64}}
    gridx = convert_length_array_units(plot_conversions, length_units, copy(gridx))
    gridy = convert_length_array_units(plot_conversions, length_units, copy(gridy))
    return gridx, gridy
end

function convert_length_array_units(
    plot_conversions::PlotConversions,
    start_units::String,
    ebarray::Vector{Float64}
)::Vector{Float64}
    conversion_func_length = get_conversion_func_length(plot_conversions, start_units)
    ebarray = apply_length_conversion_func(conversion_func_length, ebarray)
    return ebarray
end

function convert_length_units(
    plot_conversions::PlotConversions,
    start_units::String,
    length_value::Float64
)::Float64
    conversion_func_length = get_conversion_func_length(plot_conversions, start_units)
    conversion_factor = conversion_func_length(1.0)
    return length_value * conversion_factor
end

function get_conversion_func_length(
    plot_conversions::PlotConversions,
    units_start::String
)::Function
    units_end = plot_conversions.plot_units.length_units
    conv_func, units_start, units_end = get_conversion_func(
        units_start, units_end, plot_conversions.unit_conversion_data, "length_units")
    return conv_func
end

function apply_length_conversion_func(
    conv_func::Function,
    ebarray::Vector{Float64}
)::Vector{Float64}
    conversion_factor = conv_func(1.0)
    return ebarray .* conversion_factor
end

function convert_velocity_array_units(
    plot_conversions::PlotConversions,
    start_units::String,
    velocity_array::Matrix{Float64}
)::Matrix{Float64}
    conversion_func = get_conversion_func_velocity(plot_conversions, start_units)
    velocity_array = apply_conversion_factor(conversion_func, velocity_array)
    return velocity_array
end

function get_conversion_func_velocity(
    conv::PlotConversions,
    units_start::String
)::Function
    units_end = conv.plot_units.velocity_units
    conv_func, units_start, units_end = get_conversion_func(
        units_start, units_end, conv.unit_conversion_data, "velocity")
    return conv_func
end

function convert_temperature_array_units(
    plot_conversions::PlotConversions,
    units_start::String,
    ebarray2d::Matrix{Float64}
)::Matrix{Float64}
    units_end = plot_conversions.plot_units.temperature_units
    plot_conversions.plot_units.active_units = units_end
    return perform_matrix_conversion(
        plot_conversions.unit_conversion_data, units_start, units_end, ebarray2d)
end

function convert_viscosity_array_units(
    plot_conversions::PlotConversions,
    units_start::String,
    ebarray2d::Matrix{Float64}
)::Matrix{Float64}
    units_end = plot_conversions.plot_units.viscosity_units
    plot_conversions.plot_units.active_units = units_end
    return perform_matrix_conversion(
        plot_conversions.unit_conversion_data, units_start, units_end, ebarray2d)
end

function convert_strainrate_array_units(
    plot_conversions::PlotConversions,
    units_start::String,
    ebarray2d::Matrix{Float64}
)::Matrix{Float64}
    units_end = plot_conversions.plot_units.strainrate_units
    plot_conversions.plot_units.active_units = units_end
    return perform_matrix_conversion(
        plot_conversions.unit_conversion_data, units_start, units_end, ebarray2d)
end

function convert_stress_array_units(
    plot_conversions::PlotConversions,
    units_start::String,
    ebarray2d::Matrix{Float64}
)::Matrix{Float64}
    units_end = plot_conversions.plot_units.stress_units
    plot_conversions.plot_units.active_units = units_end
    return perform_matrix_conversion(
        plot_conversions.unit_conversion_data, units_start, units_end, ebarray2d)
end

function convert_pressure_array_units(
    plot_conversions::PlotConversions,
    units_start::String,
    ebarray2d::Matrix{Float64}
)::Matrix{Float64}
    units_end = plot_conversions.plot_units.pressure_units
    plot_conversions.plot_units.active_units = units_end
    return perform_matrix_conversion(
        plot_conversions.unit_conversion_data, units_start, units_end, ebarray2d)
end

function perform_matrix_conversion(
    unit_conversion_data::UnitConversionData,
    units_start::String,
    units_end::String,
    ebarray2d::Matrix{Float64}
)::Matrix{Float64}
    conversion_func = get_conversion_func_scalar(unit_conversion_data, units_start, units_end)
    ebarray2d = apply_conversion_func(conversion_func, ebarray2d)
    return ebarray2d
end

function get_conversion_func_scalar(
    unit_conversion_data::UnitConversionData,
    units_start::String,
    units_end::String
)::Function
    conv_func, units_start, units_end = get_conversion_func(
        units_start, units_end, unit_conversion_data, "scalar_units")
    return conv_func
end

function apply_conversion_factor(
    conv_func::Function,
    ebarray1d::Vector{Float64}
)::Vector{Float64}
    conversion_factor = conv_func(1.0)
    return ebarray1d .* conversion_factor
end

function apply_conversion_factor(
    conv_func::Function,
    ebarray2d::Matrix{Float64}
)::Matrix{Float64}
    conversion_factor = conv_func(1.0)
    return ebarray2d .* conversion_factor
end

function apply_conversion_func(
    conv_func::Function,
    ebarray1d::Vector{Float64}
)::Vector{Float64}
    return map(conv_func, ebarray1d)
end

function apply_conversion_func(
    conv_func::Function,
    ebarray2d::Matrix{Float64}
)::Matrix{Float64}
    return map(conv_func, ebarray2d)
end

end # module