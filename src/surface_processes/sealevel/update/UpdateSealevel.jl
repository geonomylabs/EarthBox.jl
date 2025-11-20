module UpdateSealevel

import EarthBox.PrintFuncs: print_info
import EarthBox: ConversionFuncs
import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit
import ..Options: option_names
import ..RelativeBaseLevel
import ..RelativeBaseLevel: ReferenceLithosphere
import ..RelativeBaseLevel: DensityProps
import ..RelativeBaseLevel: LithostaticPressure

function update_sealevel!(model::ModelData, ::Val{option_names.AveragePressure})::Nothing
    @timeit_memit "Finished updating sealevel using average pressure and reference lithosphere" begin
        update_sealevel_average_pressure!(model)
    end
    return nothing
end

function update_sealevel!(model::ModelData, ::Val{option_names.LeftEdge})::Nothing
    @timeit_memit "Finished updating sealevel using left edge" begin
        update_sealevel_left_edge!(model)
    end
    return nothing
end

function update_sealevel!(model::ModelData, ::Val{option_names.Constant})::Nothing
    @timeit_memit "Finished updating sealevel using fixed sea level" begin
        update_sealevel_fixed!(model)
    end
    return nothing
end

""" Update sea level using average pressure and reference lithosphere.

Sea level is defined as follows:

    y_sealevel = ymin + relative_base_level + base_level_shift

where ymin is the y-coordinate of the top of the model domain,
relative_base_level is the base level relative to ymin, and
base_level_shift is a user specified shift that can be used to model air
filled basins and dynamic topography.
"""
function update_sealevel_average_pressure!(model::ModelData)::Nothing
    ymin = model.grids.parameters.geometry.ymin.value

    relative_base_level = 
        calculate_relative_base_level_using_average_pressure(model)

    base_level_shift_time_dependent = get_time_dependent_base_level_shift(model)

    y_sealevel = ymin + relative_base_level + base_level_shift_time_dependent

    print_info("Updated sealevel (m) : $y_sealevel", level=2)
    model.topography.parameters.sealevel.y_sealevel.value = y_sealevel
    return nothing
end

function get_time_dependent_base_level_shift(model::ModelData)::Float64
    sealevel = model.topography.parameters.sealevel
    base_level_shift_meters = sealevel.base_level_shift.value
    base_level_shift_end_time_myr = sealevel.base_level_shift_end_time.value

    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_yr = ConversionFuncs.seconds_to_years(timesum)
    timesum_myr = timesum_yr/1.0e6

    base_level_shift_time_dependent = 0.0
    if timesum_myr < base_level_shift_end_time_myr
        base_level_shift_time_dependent = base_level_shift_meters
    end

    return base_level_shift_time_dependent
end

""" Calculate relative base level using average basal pressure.

The relative base level is defined with respect to the top of the model domain.

For cases where sea level is located above the sticky-rock interface,
the relative base level is negative. For cases where sea level is located
below the sticky-rock interface, the relative base level is positive.

The y-location of sea level is related to the relative base level by:

    y_sealevel = y_topo_left_edge + relative_base_level
"""
function calculate_relative_base_level_using_average_pressure(
    model::ModelData;
)::Float64
    ysize = model.grids.parameters.geometry.ysize.value

    ref_lith = model.topography.parameters.reference_lithosphere
    
    reference_lith_thickness = ReferenceLithosphere.LithosphereThicknesses(
        ysize,
        ref_lith.thickness_upper_continental_crust_ref.value,
        ref_lith.thickness_lower_continental_crust_ref.value,
        ref_lith.thickness_lithosphere_ref.value,
        ref_lith.gridy_spacing_ref.value
    )

    steady_state = model.heat_equation.parameters.steady_state
    reference_lith_thermal = ReferenceLithosphere.LithosphereThermalProps(
        ref_lith.temperature_top_ref.value,
        ref_lith.temperature_moho_ref.value,
        ref_lith.temperature_base_lith_ref.value,
        ref_lith.adiabatic_gradient_ref.value,
        steady_state.conductivity_upper_crust.value,
        steady_state.conductivity_lower_crust.value,
        steady_state.conductivity_mantle.value,
        steady_state.heat_production_upper_crust.value,
        steady_state.heat_production_lower_crust.value,
        steady_state.heat_production_mantle.value,
        steady_state.thick_thermal_lithosphere.value
    )

    gridx_pr = model.grids.arrays.pressure.gridx_pr.array
    gridy_pr = model.grids.arrays.pressure.gridy_pr.array
    pressure_array = model.stokes_continuity.arrays.pressure.pr1.array

    xnum_pr = size(gridx_pr, 1)
    ynum_pr = size(gridy_pr, 1)
    ymax_pressure = gridy_pr[ynum_pr]

    model_column_thickness_meters = ymax_pressure

    density_model = DensityProps.get_density_model(model)
    density_air = DensityProps.get_density_sticky_air(model)

    pressure_at_base_of_model_column = 0.0
    pressure_min = 1e39
    pressure_max = -1e39

    for j in 1:xnum_pr
        for i in 1:ynum_pr
            pressure = pressure_array[i, j]
            if i == ynum_pr
                pressure_at_base_of_model_column += pressure
                if pressure < pressure_min
                    pressure_min = pressure
                end
                if pressure > pressure_max
                    pressure_max = pressure
                end
            end
        end
    end

    pressure_at_base_of_model_column /= xnum_pr
    gravity_y = model.gravity.parameters.gravity_y.value
    iuse_linear_segments = ref_lith.iuse_linear_segments.value

    relative_base_level, _, _, _, _ = 
        RelativeBaseLevel.calculate_isostatic_relative_base_level_average_pressure(
            model_column_thickness_meters,
            pressure_at_base_of_model_column,
            reference_lith_thickness,
            reference_lith_thermal,
            density_model,
            density_air          = density_air,
            gravity_m_s_s        = gravity_y,
            iuse_linear_segments = iuse_linear_segments
        )

    return relative_base_level
end

""" Update sea level using left edge and reference lithosphere.

Sea level is defined as follows:

    y_sealevel = y_topo_left_edge + relative_base_level + base_level_shift

where y_topo_left_edge is the y-coordinate of the left edge of the
topography, relative_base_level is the base level relative to the
sticky-rock interface, and base_level_shift is a user specified shift that
can be used to model air filled basins and dynamic topography.
"""
function update_sealevel_left_edge!(model::ModelData)::Nothing
    sealevel = model.topography.parameters.sealevel
    base_level_shift_meters = sealevel.base_level_shift.value
    base_level_shift_end_time_myr = sealevel.base_level_shift_end_time.value

    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_yr = ConversionFuncs.seconds_to_years(timesum)
    timesum_myr = timesum_yr/1.0e6

    gridt = model.topography.arrays.gridt.array
    y_topo_left_edge = gridt[2, 1]

    relative_base_level = calculate_relative_base_level(
        model, y_topo_left_edge, iuse_linear_segments=1)

    print_info("Relative base level (m): $relative_base_level", level=2)
    print_info("y_topo_left_edge (m): $y_topo_left_edge", level=2)

    if timesum_myr < base_level_shift_end_time_myr
        y_sealevel = y_topo_left_edge + relative_base_level + base_level_shift_meters
    else
        y_sealevel = y_topo_left_edge + relative_base_level
    end

    print_info("Updated sealevel (m) : $y_sealevel", level=2)
    sealevel.y_sealevel.value = y_sealevel
    return nothing
end

""" Calculate relative base level.

The relative base level is defined with respect to sticky-rock interface.

For cases where sea level is located above the sticky-rock interface,
the relative base level is negative. For cases where sea level is located
below the sticky-rock interface, the relative base level is positive.

The y-location of sea level is related to the relative base level by:

    y_sealevel = y_topo_left_edge + relative_base_level
"""
function calculate_relative_base_level(
    model::ModelData,
    y_topo_left_edge::Float64;
    iuse_linear_segments::Int=1
)::Float64
    ysize = model.grids.parameters.geometry.ysize.value
    ref_lith = model.topography.parameters.reference_lithosphere

    reference_lith_thickness = ReferenceLithosphere.LithosphereThicknesses(
        ysize,
        ref_lith.thickness_upper_continental_crust_ref.value,
        ref_lith.thickness_lower_continental_crust_ref.value,
        ref_lith.thickness_lithosphere_ref.value,
        ref_lith.gridy_spacing_ref.value
    )

    steady_state = model.heat_equation.parameters.steady_state
    reference_lith_thermal = ReferenceLithosphere.LithosphereThermalProps(
        ref_lith.temperature_top_ref.value,
        ref_lith.temperature_moho_ref.value,
        ref_lith.temperature_base_lith_ref.value,
        ref_lith.adiabatic_gradient_ref.value,
        steady_state.conductivity_upper_crust.value,
        steady_state.conductivity_lower_crust.value,
        steady_state.conductivity_mantle.value,
        steady_state.heat_production_upper_crust.value,
        steady_state.heat_production_lower_crust.value,
        steady_state.heat_production_mantle.value,
        steady_state.thick_thermal_lithosphere.value
    )

    total_column_thickness_meters = ysize - y_topo_left_edge
    density_model = DensityProps.get_density_model(model)

    pressure_at_base_of_model_column = calculate_lithostatic_pressure_at_base_of_left_edge(
        model, y_topo_left_edge)

    relative_base_level, _, _, _, _ = 
        RelativeBaseLevel.calculate_isostatic_relative_base_level(
            total_column_thickness_meters,
            pressure_at_base_of_model_column,
            reference_lith_thickness,
            reference_lith_thermal,
            density_model,
            iuse_linear_segments=iuse_linear_segments
        )

    return relative_base_level
end

function calculate_lithostatic_pressure_at_base_of_left_edge(
    model::ModelData,
    y_topo_left_edge::Float64
)::Float64
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    marker_rho = model.markers.arrays.material.marker_rho.array

    ysize = model.grids.parameters.geometry.ysize.value

    mxstep = model.markers.parameters.distribution.mxstep.value
    mystep = model.markers.parameters.distribution.mystep.value

    _, _, pressure_gridy_from_markers = 
        LithostaticPressure.calculate_lithostatic_pressure_from_marker_swarm(
            marker_x,
            marker_y,
            marker_rho,
            ysize,
            mxstep*8,
            y_topo_left_edge, # Only consider sub-rock material
            mxstep*8,
            mystep*8
        )

    pressure_at_base_of_model_column = pressure_gridy_from_markers[end]
    return pressure_at_base_of_model_column
end

function update_sealevel_fixed!(model::ModelData)::Nothing
    sealevel = model.topography.parameters.sealevel
    base_level_shift_meters = sealevel.base_level_shift.value
    base_level_shift_end_time_myr = sealevel.base_level_shift_end_time.value

    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_yr = ConversionFuncs.seconds_to_years(timesum)
    timesum_myr = timesum_yr/1.0e6

    y_water_ini = sealevel.y_water_ini.value
    if timesum_myr < base_level_shift_end_time_myr
        shift = base_level_shift_meters
    else
        shift = 0.0
    end

    y_sealevel = y_water_ini + shift
    sealevel.y_sealevel.value = y_sealevel

    print_info = false
    if print_info
        print_sealevel_info(
            base_level_shift_end_time_myr, y_water_ini, shift, y_sealevel)
    end
    return nothing
end

function print_sealevel_info(
    base_level_shift_end_time_myr::Float64,
    y_water_ini::Float64,
    shift::Float64,
    y_sealevel::Float64
)::Nothing
    print_info("base_level_shift_end_time_myr: $base_level_shift_end_time_myr", level=2)
    print_info("y_water_ini (m): $y_water_ini", level=2)
    print_info("base_level_shift (m): $shift", level=2)
    print_info("Updated y_sealevel (m): $y_sealevel", level=2)
    return nothing
end

end # module 