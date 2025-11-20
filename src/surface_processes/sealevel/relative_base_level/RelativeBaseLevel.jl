module RelativeBaseLevel

include("core/density/DensityProps.jl")
include("core/LithostaticPressure.jl")
include("core/ReferenceLithosphere.jl")

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.PrintFuncs: print_info
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import .ReferenceLithosphere
import .DensityProps

export initialize!

const PDATA = get_eb_parameters()

struct ValidInputNames
    thickness_upper_continental_crust_ref::Symbol
    thickness_lower_continental_crust_ref::Symbol
    thickness_lithosphere_ref::Symbol
    gridy_spacing_ref::Symbol
    temperature_top_ref::Symbol
    temperature_moho_ref::Symbol
    temperature_base_lith_ref::Symbol
    adiabatic_gradient_ref::Symbol
    iuse_linear_segments::Symbol
end

"""
    initialize!(model::ModelData; kwargs...)::Nothing

Initialize reference lithosphere parameters for relative base level model. The
reference lithosphere is used to define global sea level (reference base level) as
the top of the lithosphere column in isostatic equilibrium with the average pressure
at the base of the model column [kneller23](@cite).

# Arguments
- `model::ModelData`: The model data container containing the model parameters and arrays.

# Keyword Arguments
- `$(PDATA.thickness_upper_continental_crust_ref.name)::Float64`
    - $(PDATA.thickness_upper_continental_crust_ref.description)
- `$(PDATA.thickness_lower_continental_crust_ref.name)::Float64`
    - $(PDATA.thickness_lower_continental_crust_ref.description)
- `$(PDATA.thickness_lithosphere_ref.name)::Float64`
    - $(PDATA.thickness_lithosphere_ref.description)
- `$(PDATA.gridy_spacing_ref.name)::Float64`
    - $(PDATA.gridy_spacing_ref.description)
- `$(PDATA.temperature_top_ref.name)::Float64`
    - $(PDATA.temperature_top_ref.description)
- `$(PDATA.temperature_moho_ref.name)::Float64`
    - $(PDATA.temperature_moho_ref.description)
- `$(PDATA.temperature_base_lith_ref.name)::Float64`
    - $(PDATA.temperature_base_lith_ref.description)
- `$(PDATA.adiabatic_gradient_ref.name)::Float64`
    - $(PDATA.adiabatic_gradient_ref.description)
- `$(PDATA.iuse_linear_segments.name)::Int`
    - $(PDATA.iuse_linear_segments.description)

"""
function initialize!(
    model::ModelData;
    kwargs...
)::Nothing
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
    return nothing
end

""" Calculate isostatic relative base level using average pressure.

# Arguments
- `model_column_thickness_meters`: Thickness of the model column in meters. This is 
  typically set to the height of the computational mesh.
- `avg_pressure_at_base_of_model`: Average pressure at the base of the model in Pascals.
- `reference_lith_thickness`: Reference lithosphere thicknesses.
- `reference_lith_thermal`: Reference lithosphere thermal properties.
- `density_model`: Density model for reference lithosphere.
- `density_air`: Average density of air column in kg/m^3.
- `gravity_m_s_s`: Acceleration due to gravity in m/s^2.
- `iuse_linear_segments`: If 1, use a temperature profile with four linear segments. 
    If 0, use the three-layer analytical model.

# Returns
- `relative_base_level`: Location of base level relative to basement of model column in 
  meters.
- `gridy_ref`: Reference grid y coordinates
- `temp_gridy_ref`: Reference temperature grid
- `density_gridy_ref`: Reference density grid
- `pressure_gridy_ref`: Reference pressure grid
"""
function calculate_isostatic_relative_base_level_average_pressure(
    model_column_thickness_meters::Float64,
    avg_pressure_at_base_of_model::Float64,
    reference_lith_thickness::ReferenceLithosphere.LithosphereThicknesses,
    reference_lith_thermal::ReferenceLithosphere.LithosphereThermalProps,
    density_model::DensityProps.DensityModel;
    density_air::Float64 = 1.0,
    gravity_m_s_s::Float64 = 9.8,
    iuse_linear_segments::Int64 = 1
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    
    (
        gridy_ref, 
        temp_gridy_ref, 
        density_gridy_ref, 
        pressure_gridy_ref
    ) = ReferenceLithosphere.reference_lithosphere(
            density_model,
            reference_lith_thickness,
            reference_lith_thermal;
            number_of_pressure_iterations = 3,
            iuse_linear_segments          = iuse_linear_segments,
            gravity                       = gravity_m_s_s
        )

    extra_thickness = calculate_extra_thickness(
        reference_lith_thickness.total_column_thickness_meters,
        model_column_thickness_meters
    )

    density_compensating_kgm3 = density_gridy_ref[end]

    pressure_ref_column = pressure_gridy_ref[end]

    debug = false
    if debug
        print_info("Information about the relative base level calculation", level=1)
        print_info("Pressure at the base of reference column (GPa) ", pressure_ref_column/1e9)
        print_info("Average pressure at base of model (GPa) ", avg_pressure_at_base_of_model/1e9)
        print_info("Compensating density (kg/m^3) ", density_compensating_kgm3)
    end

    if pressure_ref_column > avg_pressure_at_base_of_model
        relative_base_level = (pressure_ref_column - avg_pressure_at_base_of_model) /
            (density_compensating_kgm3 * gravity_m_s_s - density_air * gravity_m_s_s)
        relative_base_level = relative_base_level - extra_thickness
    else

        print_info("Information about the relative base level calculation", level=1)

        print_info("Ref lith array info", level=1)
        ilast = length(gridy_ref)
        print_info("Grid y coordinates (m): min = $(minimum(gridy_ref)), max = $(maximum(gridy_ref)), last = $(gridy_ref[ilast])", level=2)
        print_info("Temperature (K): min = $(minimum(temp_gridy_ref)), max = $(maximum(temp_gridy_ref)), last = $(temp_gridy_ref[end])", level=2)
        print_info("Density (kg/m^3): min = $(minimum(density_gridy_ref)), max = $(maximum(density_gridy_ref)), last = $(density_gridy_ref[ilast])", level=2)
        print_info("Pressure (Pa): min = $(minimum(pressure_gridy_ref)), max = $(maximum(pressure_gridy_ref)), last = $(pressure_gridy_ref[ilast])", level=2)

        print_info("Standard Density", level=1)
        print_info("upr_cont_crust standard density $(density_model.upr_cont_crust.standard_density)", level=2)
        print_info("lwr_cont_crust standard density $(density_model.lwr_cont_crust.standard_density)", level=2)
        print_info("upr_mantle_lithosphere standard density $(density_model.upr_mantle_lithosphere.standard_density)", level=2)
        print_info("mid_mantle_lithosphere standard density $(density_model.mid_mantle_lithosphere.standard_density)", level=2)
        print_info("lwr_mantle_lithosphere standard density $(density_model.lwr_mantle_lithosphere.standard_density)", level=2)
        print_info("asthenosphere standard density $(density_model.asthenosphere.standard_density)", level=2)

        print_info("Expansivity", level=1)
        print_info("upr_cont_crust expansivity $(density_model.upr_cont_crust.expansivity)", level=2)
        print_info("lwr_cont_crust expansivity $(density_model.lwr_cont_crust.expansivity)", level=2)
        print_info("upr_mantle_lithosphere expansivity $(density_model.upr_mantle_lithosphere.expansivity)", level=2)
        print_info("mid_mantle_lithosphere expansivity $(density_model.mid_mantle_lithosphere.expansivity)", level=2)
        print_info("lwr_mantle_lithosphere expansivity $(density_model.lwr_mantle_lithosphere.expansivity)", level=2)
        print_info("asthenosphere expansivity $(density_model.asthenosphere.expansivity)", level=2)

        print_info("Compressibility", level=1)
        print_info("upr_cont_crust compressibility $(density_model.upr_cont_crust.compressibility)", level=2)
        print_info("lwr_cont_crust compressibility $(density_model.lwr_cont_crust.compressibility)", level=2)
        print_info("upr_mantle_lithosphere compressibility $(density_model.upr_mantle_lithosphere.compressibility)", level=2)
        print_info("mid_mantle_lithosphere compressibility $(density_model.mid_mantle_lithosphere.compressibility)", level=2)
        print_info("lwr_mantle_lithosphere compressibility $(density_model.lwr_mantle_lithosphere.compressibility)", level=2)
        print_info("asthenosphere compressibility $(density_model.asthenosphere.compressibility)", level=2)

        print_info("Reference Lithosphere Thicknesses", level=1)
        print_info("total_column_thickness_meters $(reference_lith_thickness.total_column_thickness_meters)", level=2)
        print_info("thickness_upr_cont_crust_meters $(reference_lith_thickness.thickness_upr_cont_crust_meters)", level=2)
        print_info("thickness_lwr_cont_crust_meters $(reference_lith_thickness.thickness_lwr_cont_crust_meters)", level=2)
        print_info("thickness_lithosphere_meters $(reference_lith_thickness.thickness_lithosphere_meters)", level=2)
        print_info("dy_meters $(reference_lith_thickness.dy_meters)", level=2)

        print_info("Reference Lithosphere Thermal Properties", level=1)
        print_info("temperature_top_celsius $(reference_lith_thermal.temperature_top_celsius)", level=2)
        print_info("temperature_moho_celsius $(reference_lith_thermal.temperature_moho_celsius)", level=2)
        print_info("temperature_base_lith_celsius $(reference_lith_thermal.temperature_base_lith_celsius)", level=2)
        print_info("adiabatic_gradient_kelvin_km $(reference_lith_thermal.adiabatic_gradient_kelvin_km)", level=2)
        print_info("conductivity_upper_crust $(reference_lith_thermal.conductivity_upper_crust)", level=2)
        print_info("conductivity_lower_crust $(reference_lith_thermal.conductivity_lower_crust)", level=2)
        print_info("conductivity_mantle $(reference_lith_thermal.conductivity_mantle)", level=2)
        print_info("heat_production_upper_crust $(reference_lith_thermal.heat_production_upper_crust)", level=2)
        print_info("heat_production_lower_crust $(reference_lith_thermal.heat_production_lower_crust)", level=2)
        print_info("heat_production_mantle $(reference_lith_thermal.heat_production_mantle)", level=2)
        print_info("thickness_thermal_lithosphere $(reference_lith_thermal.thickness_thermal_lithosphere)", level=2)

        print_info("Other Parameters", level=1)
        print_info("model column thickness (m) $model_column_thickness_meters", level=2)
        print_info("reference column thickness (m) $(reference_lith_thickness.total_column_thickness_meters)", level=2)
        print_info("Pressure at the base of reference column (GPa) $(pressure_ref_column/1e9)", level=2)
        print_info("Average pressure at base of model (GPa) $(avg_pressure_at_base_of_model/1e9)", level=2)
        print_info("Compensating density (kg/m^3) $density_compensating_kgm3", level=2)
        print_info("Extra thickness (m) $extra_thickness", level=2)

        print_info("Ref Lith Output Arrays", level=1)
        print_info("Grid y coordinates (m): min = $(minimum(gridy_ref)), max = $(maximum(gridy_ref))", level=2)
        print_info("Temperature (K): min = $(minimum(temp_gridy_ref)), max = $(maximum(temp_gridy_ref))", level=2)
        print_info("Density (kg/m^3): min = $(minimum(density_gridy_ref)), max = $(maximum(density_gridy_ref))", level=2)
        print_info("Pressure (Pa): min = $(minimum(pressure_gridy_ref)), max = $(maximum(pressure_gridy_ref))", level=2)

        error("The average pressure at the base of the model is greater than the " *
              "pressure at the base of the reference lithosphere. Something is " *
              "wrong since the reference lithosphere has a thickness equal to " *
              "model height.")
    end

    return (relative_base_level, gridy_ref, temp_gridy_ref, 
            density_gridy_ref, pressure_gridy_ref)
end

""" Calculate isostatic relative base level.

# Arguments
- `model_column_thickness_meters`: Thickness of the model column in meters. This is 
  typically set to the height of the computational mesh.
- `pressure_at_base_of_model_column_pascals`: Pressure at the base of the model column 
  in Pascals.
- `reference_lith_thickness`: Reference lithosphere thicknesses.
- `reference_lith_thermal`: Reference lithosphere thermal properties.
- `density_model`: Density model for reference lithosphere.
- `density_water_kgm3`: Density of water in kg/m^3.
- `gravity_m_s_s`: Acceleration due to gravity in m/s^2.
- `iuse_linear_segments`: If true, use a temperature profile with four linear segments. 
  If false, use the three-layer analytical model.

# Returns
- `relative_base_level`: Location of base level relative to basement of model column in 
  meters.
- `gridy_ref`: Reference grid y coordinates
- `temp_gridy_ref`: Reference temperature grid
- `density_gridy_ref`: Reference density grid
- `pressure_gridy_ref`: Reference pressure grid
"""
function calculate_isostatic_relative_base_level(
    model_column_thickness_meters::Float64,
    pressure_at_base_of_model_column_pascals::Float64,
    reference_lith_thickness::ReferenceLithosphere.LithosphereThicknesses,
    reference_lith_thermal::ReferenceLithosphere.LithosphereThermalProps,
    density_model::DensityProps.DensityModel;
    density_water_kgm3::Float64 = 1000.0,
    gravity_m_s_s::Float64 = 9.8,
    iuse_linear_segments::Int64 = 1
)::Tuple{Float64, Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    
    (
        gridy_ref, temp_gridy_ref, 
        density_gridy_ref, pressure_gridy_ref
    ) = ReferenceLithosphere.reference_lithosphere(
            density_model,
            reference_lith_thickness,
            reference_lith_thermal;
            number_of_pressure_iterations = 3,
            iuse_linear_segments = iuse_linear_segments
        )

    extra_thickness = calculate_extra_thickness(
        reference_lith_thickness.total_column_thickness_meters,
        model_column_thickness_meters
    )

    density_compensating_kgm3 = density_gridy_ref[end]

    delta_pressure = extra_thickness * density_compensating_kgm3 * gravity_m_s_s

    pressure_model_column_update = pressure_at_base_of_model_column_pascals + delta_pressure

    pressure_ref_column = pressure_gridy_ref[end]

    if pressure_ref_column < pressure_model_column_update
        relative_base_level = (pressure_model_column_update - pressure_ref_column) /
            (density_compensating_kgm3 * gravity_m_s_s - density_water_kgm3 * gravity_m_s_s)
        # Negative value because the base level is above the top of the model
        # column and y increases with depth
        relative_base_level = -relative_base_level
    else
        # Here base level will be below the top of the model column
        relative_base_level = (pressure_ref_column - pressure_model_column_update) /
            (density_compensating_kgm3 * gravity_m_s_s)
    end

    return (relative_base_level, gridy_ref, temp_gridy_ref, 
            density_gridy_ref, pressure_gridy_ref)
end

""" Calculate extra thickness of the reference column relative to the model.

# Arguments
- `reference_column_thickness_meters`: Thickness of reference column in meters
- `model_column_thickness_meters`: Thickness of model column in meters

# Returns
- `Float64`: Extra thickness in meters
"""
function calculate_extra_thickness(
    reference_column_thickness_meters::Float64,
    model_column_thickness_meters::Float64
)::Float64
    extra_thickness = reference_column_thickness_meters - model_column_thickness_meters
    extra_thickness = max(0.0, extra_thickness)
    return extra_thickness
end

end # module