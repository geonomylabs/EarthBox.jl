"""
    ParameterRegistry

Global parameter registry containing all parameters from all model subsystems.

"""
module ParameterRegistry

import EarthBox.TupleTools: merge_named_tuples
import EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
import EarthBox.EarthBoxDtypes: InputParamDictType
import EarthBox.EarthBoxDtypes: InputMaterialParamDictType

function get_eb_parameters()::NamedTuple
    return merge_named_tuples(
        get_benchmark_parameters(),
        get_carbonate_parameters(),
        get_conversion_parameters(),
        get_gravity_parameters(),
        get_boundary_condition_parameters(),
        get_geometry_parameters(),
        get_grid_parameters(),
        get_heat_equation_parameters(),
        get_interpolation_parameters(),
        get_markers_parameters(),
        get_materials_parameters(),
        get_melting_parameters(),
        get_stokes_continuity_parameters(),
        get_timestep_parameters(),
        get_topography_parameters(),
        get_material_model_parameters(),
        get_material_override_parameters(),
        get_other_parameters()
    )
end

""" Make a dictionary where keys are parameter names and values are Parameter objects.
"""
function get_parameter_dict()::Dict{String, Union{ParameterFloat, ParameterInt, ParameterStr}}
    param_dict = Dict{String, Union{ParameterFloat, ParameterInt, ParameterStr}}()
    eb_parameters = get_eb_parameters()
    for (field_name, param) in pairs(eb_parameters)
        param_dict[param.name] = param
    end
    return param_dict
end

"""
    get_input_parameter_name_list(),

Return a tuple containing:
1. A vector of parameter object names
2. A dictionary mapping parameter names to their type information

Parameter names from input files are checked against this list.

Each dictionary key is a model parameter. Each value is a tuple containing:
(conversion_function, default_value),

where conversion_function is a Julia function that converts string values to the 
appropriate type and the default value is the value used if input is not provided.

Note that the second element of the tuple is the default value for the
parameter. However, this second element is no longer used.
"""
function get_input_parameter_name_list()::Tuple{Vector{String}, InputParamDictType}
    # Create dictionary by mapping each field to a tuple of (parser, value, units)
    eb_parameters = get_eb_parameters()
    input_param_dict = InputParamDictType(p.name => (p.parser, p.value, p.units) for (_, p) in pairs(eb_parameters))
    return collect(keys(input_param_dict)), input_param_dict
end

"""
    check_against_master(param_name::String),

Check input parameter name against master list.
Returns the conversion function for the parameter if it exists.
Throws an error if the parameter name is not valid.
"""
function check_against_master(param_name::String)::Function
    input_param_names, input_param_dict = get_input_parameter_name_list()
    if !(param_name in input_param_names)
        throw(ArgumentError("$param_name is not a valid input parameter name"))
    end
    conv_func = input_param_dict[param_name][1]
    return conv_func
end

""" Get list of missing parameter names that are in master list but not in param_dict and not obsolete.
"""
function get_missing_names(
    param_dict::InputParamDictType,
    master_param_names::Vector{String},
    obsolete_param_names::Vector{String}
)
    missing_names = String[]
    for param_name in master_param_names
        if !(param_name in keys(param_dict)) && !(param_name in obsolete_param_names)
            push!(missing_names, param_name)
        end
    end
    return missing_names
end

function get_master_material_parameters()::Tuple{Vector{String}, InputMaterialParamDictType}
    eb_parameters = get_material_model_parameters()
    master_parameters = InputMaterialParamDictType(
        p.name => (p.parser, p.value, p.units) for (_, p) in pairs(eb_parameters)
    )
    master_names = collect(keys(master_parameters))
    return master_names, master_parameters
end

"""
    check_against_material_master(param_name::String)

Check input parameter name against material master list.
Returns the conversion function for the parameter if it exists.
Throws an error if the parameter name is not valid.
"""
function check_against_material_master(param_name::String)::Function
    master_names, master_parameters = get_master_material_parameters()
    if !(param_name in master_names)
        throw(ArgumentError("$param_name is not a valid material parameter name"))
    end
    conv_func = master_parameters[param_name][1]
    return conv_func
end

""" Get list of missing parameter names that are in master list but not in param_dict and not obsolete.
"""
function get_missing_material_names(
    param_dict::InputMaterialParamDictType,
    master_param_names::Vector{String},
    obsolete_param_names::Vector{String}
)
    missing_names = String[]
    for param_name in master_param_names
        if !(param_name in keys(param_dict)) && !(param_name in obsolete_param_names)
            push!(missing_names, param_name)
        end
    end
    return missing_names
end

function get_benchmark_parameters()::NamedTuple
    return (
        iuse_ramberg_post_processing = ParameterInt(
            0, "iuse_ramberg_post_processing", "None", 
            "Calculate quantities for Rayleigh-Taylor benchmark: 0 = off; 1 = on"),
        iuse_viscous_block_processing = ParameterInt(
            0, "iuse_viscous_block_processing", "None", 
            "Calculate quantities for viscous block benchmark: 0 = off; 1 = on"),
        iuse_conbox_post_processing = ParameterInt(
            0, "iuse_conbox_post_processing", "None", 
            "Post-processing for convection in a box: 0 = off; 1 = on"),
    )
end

function get_carbonate_parameters()::NamedTuple
    return (
        iuse_carb = ParameterInt(
            0, "iuse_carb", "None", "0 = uniform deposition; 1 = deposition on highs"),
        photic_thick_m = ParameterFloat(
            0.0, "photic_thick_m", "meters", "Depth of photic zone"),
        carb_growth_rad = ParameterFloat(
            0.0, "carb_growth_rad", "meters", "Carbonate growth radius"),
        carb_base_rate = ParameterFloat(
            0.0, "carb_base_rate", "fraction", "Probability of nucleation"),
        carb_time_myr = ParameterFloat(
            0.0, "carb_time_myr", "Myr", "Carbonate growth delay time relative to model start"),
        carb_jump_time_myr = ParameterFloat(
            0.0, "carb_jump_time_myr", "Myr", "Time of jump in carbonate base rate"),
        carb_base_rate_jump = ParameterFloat(
            0.0, "carb_base_rate_jump", "fraction", "Carbonate base rate jump"),
    )
end

function get_conversion_parameters()::NamedTuple
    return (
        KtoC = ParameterFloat(
            273.0, "KtoC", "None", "Conversion factor for Kelvins to Celsius: T_C = T_K - KtoC"),
        sec_per_yr = ParameterFloat(
            365.25*24.0*3600.0, "sec_per_yr", "s/yr", "Seconds per year"),
        sec_per_Myr = ParameterFloat(
            365.25*24.0*3600.0*1e+6, "sec_per_Myr", "s/Myr", "Seconds per million years"),
        cm_yr2m_s = ParameterFloat(
            1.0/(100.0*365.25*24.0*3600.0), "cm_yr2m_s", "m/s/cm/yr", 
            "Conversion factor for cm/yr to m/s"
            ),
    )
end

function get_gravity_parameters()::NamedTuple
    return (
        gravity_x = ParameterFloat(
            0.0, "gravity_x", "m/s/s", "x-component of gravity acceleration in m/s/s"),
        gravity_y = ParameterFloat(
            9.8, "gravity_y", "m/s/s", "y-component of gravity acceleration in m/s/s"),
        turn_off_gravity_y = ParameterInt(
            0, "turn_off_gravity_y", "None", "Flag to turn off y-component of gravity: 0 = off; 1 = on"),
        nsteps_turn_off_gravity = ParameterInt(
            100, "nsteps_turn_off_gravity", "None", "Number of steps to gradually turn off gravity"),
    )
end

function get_boundary_condition_parameters()::NamedTuple
    return (
        # BCs - Pressure
        pressure_bc = ParameterFloat(
            10000.0, "pressure_bc", "Pa", "Pressure boundary condition in Pa."),

        # BCs - Velocity parameters
        velocity = ParameterFloat(
            0.0, "velocity", "m/s", "General velocity boundary condition"),
        full_velocity_extension = ParameterFloat(
            0.0, "full_velocity_extension", "m/s", "Full extension rate in m/s"),
        full_velocity_extension_step1 = ParameterFloat(
            NaN, "full_velocity_extension_step1", "m/s", "Full extension rate at step 1 in m/s"),
        full_velocity_extension_step2 = ParameterFloat(
            NaN, "full_velocity_extension_step2", "m/s", "Full extension rate at step 2 in m/s"),
        full_velocity_contraction = ParameterFloat(
            0.0, "full_velocity_contraction", "m/s", "Full contraction rate in m/s"),
        velocity_shear = ParameterFloat(
            0.0, "velocity_shear", "m/s", "Shear velocity in m/s"),
        velocity_rotation = ParameterFloat(
            0.0, "velocity_rotation", "m/s", "Rotation velocity in m/s"),
        iuse_strain_rate = ParameterInt(
            0, "iuse_strain_rate", "None", 
            "0 use full extension velocity to define velocity bc's; 1 use strain rate"
            ),
        strain_rate_bc = ParameterFloat(
            1e-15, "strain_rate_bc", "1/s", "Strain rate boundary condition in 1/s"),
        vyu = ParameterFloat(
            0.0, "vyu", "m/s", "Velocity in y direction at upper boundary in m/s"),
        vyl = ParameterFloat(
            0.0, "vyl", "m/s", "Velocity in y direction at lower boundary in m/s"),
        velocity_internal_x = ParameterFloat(
            0.0, "velocity_internal_x", "m/s", "X component of internal velocity in m/s"),
        velocity_internal_y = ParameterFloat(
            0.0, "velocity_internal_y", "m/s", "Y component of internal velocity in m/s"),
        plate_thickness = ParameterFloat(
            100000.0, "plate_thickness", "m",
            "Thickness of plate in meters used for inflow-outflow depth-dependent extension and "
            *"contraction boundary conditions"
            ),
        smoothing_thickness = ParameterFloat(
            10000.0, "smoothing_thickness", "m", 
            "Thickness in meters for velocity smoothing between inflow and outflow boundaries"
            ),
        velocity_inflow_left = ParameterFloat(
            0.0, "velocity_inflow_left", "m/s", "Inflow velocity at left boundary in m/s"),
        velocity_inflow_right = ParameterFloat(
            0.0, "velocity_inflow_right", "m/s", "Inflow velocity at right boundary in m/s"),
        velocity_inflow_smooth_avg_left = ParameterFloat(
            0.0, "velocity_inflow_smooth_avg_left", "m/s", 
            "Smoothed average inflow velocity at left boundary in m/s"
            ),
        velocity_inflow_smooth_avg_right = ParameterFloat(
            0.0, "velocity_inflow_smooth_avg_right", "m/s", 
            "Smoothed average inflow velocity at right boundary in m/s"
            ),
        iuse_velocity_stop = ParameterInt(
            0, "iuse_velocity_stop", "None", "Flag for velocity stopping: 0 = off; 1 = on"),
        velocity_stop_time = ParameterFloat(
            0.0, "velocity_stop_time", "Myr", "Time at which velocity stops in Myr"),
        ivelocity_stop_counter = ParameterInt(
            0, "ivelocity_stop_counter", "None", "Counter for velocity stop: 0 = off; 1 = on"),

        # BCs - Temperature parameters
        temperature_top = ParameterFloat(
            0.0, "temperature_top", "K", "Temperature along the top boundary of the model domain in Kelvin"),
        temperature_bottom = ParameterFloat(
            0.0, "temperature_bottom", "K", "Temperature at the base of the model domain in Kelvin"),
        temperature_left = ParameterFloat(
            0.0, "temperature_left", "K", "Temperature along the left boundary of the model domain in Kelvin"),
        temperature_right = ParameterFloat(
            0.0, "temperature_right", "K", "Temperature along the right boundary of the model domain in Kelvin"),
        iuse_bottom_transient = ParameterInt(
            0, "iuse_bottom_transient", "None", 
            "Flag to use transient temperature at the bottom boundary"
            ),
        temperature_bottom_transient = ParameterFloat(
            NaN, "temperature_bottom_transient", "K", 
            "transient temperature at the bottom boundary"
            ),
        temperature_bottom_original = ParameterFloat(
            -99999.0, "temperature_bottom_original", "K", 
            "Original temperature at the bottom boundary"
            ),
        start_time_bottom_transient = ParameterFloat(
            0.0, "start_time_bottom_transient", "Myr", 
            "Start time (Myr) for transient temperature at the bottom boundary"
            ),
        end_time_bottom_transient = ParameterFloat(
            0.0, "end_time_bottom_transient", "Myr", 
            "End time (Myr) for transient temperature at the bottom boundary"
            ),
        delta_temperature_transient = ParameterFloat(
            NaN, "delta_temperature_transient", "deltaK", 
            "Temperature perturbation for initial transient temperature at the bottom boundary"
            ),
        temperature_base_lith_warmer_initial = ParameterFloat(
            NaN, "temperature_base_lith_warmer_initial", "K", 
            "Temperature at the base of the lithosphere"
            ),
        temperature_bottom_warmer_initial = ParameterFloat(
            NaN, "temperature_bottom_warmer_initial", "K", 
            "Initial elevated temperature at the bottom of the model"
            ),
        temperature_bottom_cooler_final = ParameterFloat(
            NaN, "temperature_bottom_cooler_final", "K", 
            "Cooler transient temperature at the bottom of the model"
            ),

        # BCs - Velocity step parameters
        iuse_velocity_step = ParameterInt(
            0, "iuse_velocity_step", "None", 
            "Increase/decrease full extension velocity by factor at specified time"
            ),
        velocity_step_factor = ParameterFloat(
            NaN, "velocity_step_factor", "None", 
            "Factor used to increase/decrease velocity at specified time"
            ),
        timestep_adjustment_factor = ParameterFloat(
            1.0, "timestep_adjustment_factor", "None", 
            "Factor used to adjust time step when velocity is increased/decreased"
            ),
        velocity_step_time = ParameterFloat(
            NaN, "velocity_step_time", "Myr", "Time of velocity increase/decrease in Myr"),
        velocity_second_step_time = ParameterFloat(
            NaN, "velocity_second_step_time", "Myr", "Time of second velocity increase/decrease in Myr"),
        ivelocity_step_counter = ParameterInt(
            0, "ivelocity_step_counter", "None", "Integer counter for velocity step"),
        velocity_second_step_factor = ParameterFloat(
            NaN, "velocity_second_step_factor", "None", 
            "Factor used to increase/decrease velocity at specified time"
            ),
        timestep_second_adjustment_factor = ParameterFloat(
            1.0, "timestep_second_adjustment_factor", "None", 
            "Factor used to adjust time step when velocity is increased/decreased"
            ),

        # BCs - Bc option parameters
        itype_bc = ParameterInt(
            0, "itype_bc", "None", "Boundary condition option for Stokes and heat equations"),
        stype_bc = ParameterStr(
            "None", "stype_bc", "None", "Boundary condition option name"),
        pressure_bc_mode = ParameterFloat(
            0.0, "pressure_bc_mode", "None", "Pressure bc: 0 upper left cell 1 top and bottom"),
    )
end

function get_geometry_parameters()::NamedTuple
    return (
        # *****************************
        # Geometry parameters
        # *****************************

        # Geometry - Sticky air parameters
        thick_air = ParameterFloat(
            NaN, "thick_air", "m", "Thickness in meters of sticky air in meters"),
        # Geometry - Sandbox parameters
        nsand_layers = ParameterInt(
            0, "nsand_layers", "None", 
            "Number of sand layers used in sandbox extension and shortening models"
            ),
        y_sand_air_interface = ParameterFloat(
            0.0, "y_sand_air_interface", "m", 
            "Y-location in meters of interface between sand and air along the left boundary used in "
            *"sandbox extension and shortening models"
            ),
        y_top_microbeads = ParameterFloat(
            0.0, "y_top_microbeads", "m", 
            "Y-location in meters of top of microbeads used in sandbox shortening models"
            ),
        y_bottom_microbeads = ParameterFloat(
            0.0, "y_bottom_microbeads", "m", 
            "Y-location in meters of bottom of microbeads used in sandbox shortening models"
            ),
        x_left_ramp = ParameterFloat(
            0.0, "x_left_ramp", "m", 
            "X-location in meters of left-edge of ramp used in sandbox shortening models"
            ),
        x_right_ramp = ParameterFloat(
            0.0, "x_right_ramp", "m", 
            "X-location in meters of right-edge of ramp used in sandbox shortening models"
            ),
        ramp_dip_deg = ParameterFloat(
            0.0, "ramp_dip_deg", "Degrees", 
            "Dip of ramp in degrees used in sandbox shortening models used in sandbox shortening models"
            ),
        pdms_layer_width = ParameterFloat(
            0.0, "pdms_layer_width", "m", 
            "Width in meters of PDMS (polydimethylsiloxane) layer located in the central bottom part "
            *"of model used in sandbox extension models"
            ),
        pdms_layer_thickness = ParameterFloat(
            0.0, "pdms_layer_thickness", "m", 
            "Thickness in meters of PDMS (polydimethylsiloxane) layer used in sandbox extension models"
            ),

        # Geometry - Crustal hole parameters
        xhole_start = ParameterFloat(
            1000.0, "xhole_start", "m", "x-coordinate in meters of the start of the crustal hole"),
        xhole_middle = ParameterFloat(
            1100.0, "xhole_middle", "m", "x-coordinate in meters of the middle of the crustal hole"),
        xhole_end = ParameterFloat(
            1200.0, "xhole_end", "m", "x-coordinate in meters of the end of the crustal hole"),
        xhole_depth = ParameterFloat(
            10.0, "xhole_depth", "m", "Depth in meters of the crustal hole"),

        # Geometry - Mobile wall parameters
        x_left_mobile_wall = ParameterFloat(
            0.0, "x_left_mobile_wall", "m", 
            "X-location in meters of left-edge of mobile wall used in sandbox extension and shortening models"),
        x_right_mobile_wall = ParameterFloat(
            0.0, "x_right_mobile_wall", "m", 
            "X-location in meters of right-edge of mobile wall used in sandbox extension and shortening models"),
        y_top_mobile_wall = ParameterFloat(
            0.0, "y_top_mobile_wall", "m", 
            "Y-location in meters of top-edge of mobile wall used in sandbox extension and shortening models"),
        y_bottom_mobile_wall = ParameterFloat(
            0.0, "y_bottom_mobile_wall", "m", 
            "Y-location in meters of bottom-edge of mobile wall used in sandbox extension and shortening models"),
        plate_extension_width = ParameterFloat(
            0.0, "plate_extension_width", "m", 
            "Width in meters of plate extension used in sandbox extension models"),
        plate_extension_thickness = ParameterFloat(
            0.0, "plate_extension_thickness", "m", 
            "Thickness in meters of plate extension used in sandbox extension models"),

        # Geometry - Rayleigh Taylor parameters
        depth_interface_h1 = ParameterFloat(
            100000.0, "depth_interface_h1", "m", 
            "Depth in meters of the interface for Rayleigh-Taylor instability"
            ),
        wave_length_lambda = ParameterFloat(
            100000.0, "wave_length_lambda", "m", 
            "Wavelength in meters of the perturbation for Rayleigh-Taylor instability"
            ),
        amplitude_initial = ParameterFloat(
            500.0, "amplitude_initial", "m", 
            "Initial amplitude in meters of the perturbation for Rayleigh-Taylor instability"
            ),

        # Geometry - Earth layering parameters
        thick_upper_crust = ParameterFloat(
            NaN, "thick_upper_crust", "m", "Thickness in meters of upper crust"),
        thick_lower_crust = ParameterFloat(
            NaN, "thick_lower_crust", "m", "Thickness in meters of lower crust"),
        thick_upper_lith = ParameterFloat(
            NaN, "thick_upper_lith", "m", "Thickness in meters of upper mantle lithosphere"),
        thick_middle_lith = ParameterFloat(
            NaN, "thick_middle_lith", "m", "Thickness in meters of middle lithosphere"),
        thick_lower_lith = ParameterFloat(
            NaN, "thick_lower_lith", "m", "Thickness in meters of lower mantle lithosphere"),
        thick_lith = ParameterFloat(
            NaN, "thick_lith", "m", "Thickness in meters of lithosphere"),
        thick_crust = ParameterFloat(
            NaN, "thick_crust", "m", "Thickness in meters of crust"),
        thick_mantle_lith = ParameterFloat(
            NaN, "thick_mantle_lith", "m", "Thickness in meters of mantle lithosphere"),

        # Geometry - Fracture zone parameters
        sediment_thickness = ParameterFloat(
            0.0, "sediment_thickness", "m", "Thickness in meters of sediment layer"),
        basaltic_oceanic_crust_thickness = ParameterFloat(
            0.0, "basaltic_oceanic_crust_thickness", "m", 
            "Thickness in meters of basaltic oceanic crust layer"
            ),
        gabbroic_oceanic_crust_thickness = ParameterFloat(
            0.0, "gabbroic_oceanic_crust_thickness", "m", 
            "Thickness in meters of gabbroic oceanic crust layer"
            ),
        thickness_of_younger_lithosphere = ParameterFloat(
            0.0, "thickness_of_younger_lithosphere", "m", "Thickness in meters of younger lithosphere"),
        thickness_of_older_lithosphere = ParameterFloat(
            0.0, "thickness_of_older_lithosphere", "m", "Thickness in meters of older lithosphere"),
        thickness_of_weak_lithosphere = ParameterFloat(
            0.0, "thickness_of_weak_lithosphere", "m", "Thickness in meters of weak lithosphere"),
        x_fracture_zone_start = ParameterFloat(
            0.0, "x_fracture_zone_start", "m", "Start of fracture zone in x-direction in meters"),
        x_fracture_zone_end = ParameterFloat(
            0.0, "x_fracture_zone_end", "m", "End of fracture zone in x-direction in meters"),

        # Geometry - Internal velocity zone parameters
        xindex_vx_internal = ParameterInt(
            -1, "xindex_vx_internal", "None", "X-index for internal vx velocity"),
        yindex_min_vx_internal = ParameterInt(
            -1, "yindex_min_vx_internal", "None", "Minimum y-index for internal vx velocity"),
        yindex_max_vx_internal = ParameterInt(
            -1, "yindex_max_vx_internal", "None", "Maximum y-index for internal vx velocity"),
        xindex_vy_internal = ParameterInt(
            -1, "xindex_vy_internal", "None", "X-index for internal vy velocity"),
        yindex_min_vy_internal = ParameterInt(
            -1, "yindex_min_vy_internal", "None", "Minimum y-index for internal vy velocity"),
        yindex_max_vy_internal = ParameterInt(
            -1, "yindex_max_vy_internal", "None", "Maximum y-index for internal vy velocity"),
        x_vx_internal = ParameterFloat(
            0.0, "x_vx_internal", "m", "X-coordinate in meters for internal vx velocity"),
        y_min_vx_internal = ParameterFloat(
            0.0, "y_min_vx_internal", "m", "Minimum y-coordinate in meters for internal vx velocity"),
        y_max_vx_internal = ParameterFloat(
            0.0, "y_max_vx_internal", "m", "Maximum y-coordinate in meters for internal vx velocity"),
        x_vy_internal = ParameterFloat(
            0.0, "x_vy_internal", "m", "X-coordinate in meters for internal vy velocity"),
        y_min_vy_internal = ParameterFloat(
            0.0, "y_min_vy_internal", "m", "Minimum y-coordinate in meters for internal vy velocity"),
        y_max_vy_internal = ParameterFloat(
            0.0, "y_max_vy_internal", "m", "Maximum y-coordinate in meters for internal vy velocity"),

        # Geometry - Litho strong zones parameters
        x_left_strong = ParameterFloat(
            100000.0, "x_left_strong", "m", 
            "x-coordinate in meters of the left edge of the lithospheric strong zone"
            ),
        x_right_strong = ParameterFloat(
            200000.0, "x_right_strong", "m", 
            "x-coordinate in meters of the right edge of the lithospheric strong zone"
            ),
        iuse_strong_zones = ParameterInt(
            0, "iuse_strong_zones", "None", "Use strong zones (0=off, 1=on)"),

        # Geometry - Plume parameters
        plume_radius = ParameterFloat(
            200000.0, "plume_radius", "m", "Radius in meters of the mantle plume"),
        plume_center_x = ParameterFloat(
            100000.0, "plume_center_x", "m", "x-coordinate in meters of the plume center"),
        plume_center_y = ParameterFloat(
            600000.0, "plume_center_y", "m", "y-coordinate in meters of the plume center"),
        plume_head_thick = ParameterFloat(
            100000.0, "plume_head_thick", "m", "Thickness in meters of the plume head"),
        delta_temperature_plume = ParameterFloat(
            100.0, "delta_temperature_plume", "K", 
            "Temperature difference of the plume relative to background mantle"
            ),
        iuse_plume = ParameterInt(
            0, "iuse_plume", "None", "Flag for using plume (0=off, 1=on)"),
        plume_start_time = ParameterFloat(
            0.0, "plume_start_time", "Myr", "Time in Myr when plume starts"),
        iplume_state = ParameterInt(
            0, "iplume_state", "None", "State of the plume: 0 - not initiated; 1 = initiated"),

        # Geometry - Weak fault parameters
        fault_dip_degrees = ParameterFloat(
            0.0, "fault_dip_degrees", "Degrees", 
            "Dip in degrees of the weak zone used to approximate a fault"
            ),
        fault_thickness = ParameterFloat(
            0.0, "fault_thickness", "m", 
            "Thickness in meters of the weak zone used to approximate a fault"
            ),
        x_initial_fault = ParameterFloat(
            0.0, "x_initial_fault", "m", 
            "Initial x-coordinate in meters of the weak zone used to approximate a fault"
            ),
        fault_height = ParameterFloat(
            0.0, "fault_height", "m", "Height in meters of the weak zone used to approximate a fault"),
        iuse_weak_fault = ParameterInt(
            0, "iuse_weak_fault", "None", "Use weak fault (0=off, 1=on)"),

        # Geometry - Weak seed parameters
        w_seed = ParameterFloat(
            1000.0, "w_seed", "m", "Width in meters of weak seed"),
        x_seed = ParameterFloat(
            100000.0, "x_seed", "m", "x-location in meters of weak seed"),
        y_seed = ParameterFloat(
            35000.0, "y_seed", "m", "y-location in meters of weak seed"),
        iuse_weak_seed = ParameterInt(
            0, "iuse_weak_seed", "None", "Use weak seed (0=off, 1=on)"),
        )
end

function get_grid_parameters()::NamedTuple
    return (
        # *****************************
        # Grid Parameters
        # *****************************

        # Grid - Geometry parameters (requires constructor args - using defaults),
        ynum = ParameterInt(
            101, "ynum", "None", "Number of basic grid points in y-direction"),
        xnum = ParameterInt(
            101, "xnum", "None", "Number of basic grid points in x-direction"),
        ysize = ParameterFloat(
            100000.0, "ysize", "m", "Height of model in meters"),
        xsize = ParameterFloat(
            100000.0, "xsize", "m", "Width of model in meters"),
        ystpavg = ParameterFloat(
            0.0, "ystpavg", "m", "Average spacing of basic grid in y-direction in meters"),
        xstpavg = ParameterFloat(
            0.0, "xstpavg", "m", "Average spacing of basic grid in x-direction in meters"),
        ymin = ParameterFloat(
            0.0, "ymin", "m", "Minimum y-coordinate of model domain in meters"),
        ymax = ParameterFloat(
            100000.0, "ymax", "m", "Maximum y-coordinate of model domain in meters"),
        xmin = ParameterFloat(
            0.0, "xmin", "m", "Minimum x-coordinate of model domain in meters"),
        xmax = ParameterFloat(
            100000.0, "xmax", "m", "Maximum x-coordinate of model domain in meters"),
        xsize_start = ParameterFloat(
            100000.0, "xsize_start", "m", "Initial width of model in meters"),

        # Grid - Options parameters
        itype_grid = ParameterInt(
            0, "itype_grid", "None", "Grid option ID"),
        stype_grid = ParameterStr(
            "None", "stype_grid", "None", "Grid option name"),

        # Grid - Refinement parameters
        dx_highres = ParameterFloat(
            NaN, "dx_highres", "m", 
            "Constant horizontal grid resolution in high-resolution area in meters"
            ),
        dx_lowres = ParameterFloat(
            NaN, "dx_lowres", "m", 
            "Average horizontal grid resolution in low-resolution area in meters"
            ),
        xo_highres = ParameterFloat(
            NaN, "xo_highres", "m", "x-location of first node of high-resolution area in meters"),
        ixo_highres = ParameterInt(
            -9999, "ixo_highres", "None", "x-index of first node of high-resolution area"),
        xf_highres = ParameterFloat(
            NaN, "xf_highres", "m", "x-location of last node of high-resolution area in meters"),
        dy_highres = ParameterFloat(
            NaN, "dy_highres", "m", "Constant vertical grid resolution in high-resolution area in meters"),
        dy_lowres = ParameterFloat(
            NaN, "dy_lowres", "m", "Average vertical grid resolution in low-resolution area in meters"),
        yf_highres = ParameterFloat(
            NaN, "yf_highres", "m", "y-location of last node of high-resolution area in meters"),
        iuse_trench = ParameterInt(
            0, "iuse_trench", "None", "Flag to use trench location for refinement"),
        trench_location = ParameterFloat(
            NaN, "trench_location", "m", "x-location of trench in meters"),
        iuse_refinement_delay = ParameterInt(
            0, "iuse_refinement_delay", "None", "Flag to delay grid refinement"),
        refinement_time = ParameterFloat(
            NaN, "refinement_time", "Myr", "Time to start grid refinement in Myr"),
        refinement_flag = ParameterInt(
            0, "refinement_flag", "None", "Flag indicating if grid is currently refined"),
        iuse_refinement_gap = ParameterInt(
            0, "iuse_refinement_gap", "None", "Flag to temporarily disable refinement"),
        refinement_gap_start_time = ParameterFloat(
            NaN, "refinement_gap_start_time", "Myr", "Time to start refinement gap in Myr"),
        refinement_gap_end_time = ParameterFloat(
            NaN, "refinement_gap_end_time", "Myr", "Time to end refinement gap in Myr"),

        )
end

function get_heat_equation_parameters()::NamedTuple
    return (
        # *****************************
        # Heat Equation Parameters
        # *****************************
        
        # Heat equation - Temp change limit parameters
        max_temp_change = ParameterFloat(
            20.0, "max_temp_change", "K", "Maximal temperature change allowed for one timestep in Kelvin"),

        # Heat equation - Steady state parameters
        thick_thermal_lithosphere = ParameterFloat(
            NaN, "thick_thermal_lithosphere", "m", "Thickness of the thermal lithosphere in meters"),
        conductivity_upper_crust = ParameterFloat(
            3.0, "conductivity_upper_crust", "W/m/K", "Thermal conductivity of upper crust W/m/K"),
        conductivity_lower_crust = ParameterFloat(
            3.0, "conductivity_lower_crust", "W/m/K", "Thermal conductivity of lower crust W/m/K"),
        conductivity_mantle = ParameterFloat(
            3.0, "conductivity_mantle", "W/m/K", "Thermal conductivity of mantle W/m/K"),
        heat_production_upper_crust = ParameterFloat(
            1.5e-6, "heat_production_upper_crust", "W/m^3", "Radiogenic heat production in upper crust W/m^3"),
        heat_production_lower_crust = ParameterFloat(
            2.0e-7, "heat_production_lower_crust", "W/m^3", "Radiogenic heat production in lower crust W/m^3"),
        heat_production_mantle = ParameterFloat(
            0.0, "heat_production_mantle", "W/m^3", "Radiogenic heat production in mantle W/m^3"),
        amplitude_perturbation = ParameterFloat(
            0.0, "amplitude_perturbation", "m", "Amplitude of central thermal perturbation in meters"),
        width_perturbation = ParameterFloat(
            0.0, "width_perturbation", "m", "Width of central thermal perturbation in meters"),
        temperature_surface = ParameterFloat(
            0.0, "temperature_surface", "K", 
            "Initial temperature at the surface for lithospheric models with sticky air in Kelvin"
            ),
        temperature_moho = ParameterFloat(
            0.0, "temperature_moho", "K", "Initial temperature at the Moho for lithospheric models in Kelvin"),
        temperature_base_lith = ParameterFloat(
            NaN, "temperature_base_lith", "K", "Initial temperature at the base of the lithosphere in Kelvin"),
        # This may be obsolete
        temperature_base_lithosphere = ParameterFloat(
            NaN, "temperature_base_lithosphere", "K", 
            "Temperature at the base of the thermal lithosphere in Kelvin"
            ),

        # Heat equation - Initial condition parameters
        itype_temp = ParameterInt(
            0, "itype_temp", "None", "Option flag for initial temperature condition"),
        stype_temp = ParameterStr(
            "None", "stype_temp", "None", "Initial temperature condition option name"),
        temperature_uniform = ParameterFloat(
            0.0, "temperature_uniform", "K", "Uniform initial temperature value in Kelvin"),
        temperature_of_wave = ParameterFloat(
            0.0, "temperature_of_wave", "K", "Temperature of wave in Kelvin"),
        age_lithosphere_left = ParameterFloat(
            0.0, "age_lithosphere_left", "Myr", "Age of lithosphere to the left of fracture zone in Myr"),
        age_lithosphere_right = ParameterFloat(
            0.0, "age_lithosphere_right", "Myr", "Age of lithosphere to the right of fracture zone in Myr"),
        thermal_lithosphere_depth_left = ParameterFloat(
            100000.0, "thermal_lithosphere_depth_left", "m", 
            "Thermal lithosphere depth to the left of fracture zone in meters"
            ),
        thermal_lithosphere_depth_right = ParameterFloat(
            200000.0, "thermal_lithosphere_depth_right", "m", 
            "Thermal lithosphere depth to the right of fracture zone in meters"
            ),
        thermal_diffusivity = ParameterFloat(
            1e-6, "thermal_diffusivity", "m^2/s", "Thermal diffusivity in m^2/s"),
        adiabatic_gradient = ParameterFloat(
            0.4, "adiabatic_gradient", "K/km", "Adiabatic gradient in K/km"),
        
        # Heat equation - RhoCp parameters
        itype_rhocp = ParameterInt(
            0, "itype_rhocp", "None", "Integer flag for density*heat_capacity model"),
        stype_rhocp = ParameterStr(
            "None", "stype_rhocp", "None", "Density*heat_capacity model option name"),
        maximum_heat_capacity = ParameterFloat(
            1200.0, "maximum_heat_capacity", "J/kg/K", 
            "Maximum heat capacity for rock properties used with Waples temperature-dependent model in J/kg/k"
            ),

        # Heat equation - Build parameters (requires constructor args - using defaults),
        ibuild_heat = ParameterInt(
            1, "ibuild_heat", "None", "0 define full system, 1 only define non-zero elements"),
        Nheat = ParameterInt(
            10201, "Nheat", "None", "Number of rows or columns in the large NxN matrix"),
        nonzero_max_heat = ParameterInt(
            316231, "nonzero_max_heat", "None", 
            "Maximum number of non-zero values allowed for the large matrix"
            ),

        # Heat equation - Heat options parameters
        iuse_heat = ParameterInt(
            1, "iuse_heat", "None", "Solve heat equation: 0 off; 1 on"),
        iuse_shear_heating = ParameterInt(
            0, "iuse_shear_heating", "None", "Shear heating: 0 off; 1 on"),
        iuse_adiabatic_heating = ParameterInt(
            0, "iuse_adiabatic_heating", "None", "Adiabatic heating: 0 off; 1 on"),
        iuse_sticky_correction = ParameterInt(
            0, "iuse_sticky_correction", "None", 
            "Sticky temperature correction (reset to top boundary temperature): 0 off; 1 on"
            ),

        # Heat equation - Thermal cond parameters
        itype_conductivity = ParameterInt(
            0, "itype_conductivity", "None", "Integer flag for thermal conductivity model"),
        stype_conductivity = ParameterStr(
            "None", "stype_conductivity", "None", "Thermal conductivity model option name"),

        )
end

function get_interpolation_parameters()::NamedTuple
    return (
        # *****************************
        # Interpolation Parameters
        # *****************************
        
        # Interpolation parameters
        iuse_initial_temp_interp = ParameterInt(
            0, "iuse_initial_temp_interp", "None", 
            "Interpolate nodal temperatures back to markers at start: 0 off; 1 on"
            ),
        iuse_harmonic_avg_normal_viscosity = ParameterInt(
            0, "iuse_harmonic_avg_normal_viscosity", "None", 
            "Harmonic averaging of shear viscosity for normal viscosity: 0 off; 1 on"
            ),
        )
end

function get_markers_parameters()::NamedTuple
    return (
        # *****************************
        # Markers Parameters
        # *****************************
        
        # Markers - Recycling parameters
        imantle = ParameterInt(
            8, "imantle", "None", "Material ID of the mantle. This is no longer used and will be removed in the future"),
        # Markers - Distribution parameters (requires constructor args - using defaults),
        iuse_random = ParameterInt(
            0, "iuse_random", "None", "Randomize initial marker distribution: 0 off; 1 on"),
        dx_marker = ParameterFloat(
            1000.0, "dx_marker", "m", "Initial average marker spacing in horizontal direction in meters"),
        dy_marker = ParameterFloat(
            1000.0, "dy_marker", "m", "Initial average marker spacing in vertical direction in meters"),
        nmarkers_cell_x = ParameterFloat(
            3.0, "nmarkers_cell_x", "None", "Number of markers per cell in horizontal direction"),
        nmarkers_cell_y = ParameterFloat(
            3.0, "nmarkers_cell_y", "None", "Number of markers per cell in vertical direction"),
        mxnum = ParameterInt(
            300, "mxnum", "None", "Total number of markers in horizontal direction"),
        mynum = ParameterInt(
            300, "mynum", "None", "Total number of markers in vertical direction"),
        marknum = ParameterInt(
            90000, "marknum", "None", "Number of markers"),
        mxstep = ParameterFloat(
            333.33, "mxstep", "m", "Distance between markers in horizontal direction in meters"),
        mystep = ParameterFloat(
            333.33, "mystep", "m", "Distance between markers in vertical direction in meters"),

        # Markers - Advection parameters
        marker_cell_displ_max = ParameterFloat(
            0.5, "marker_cell_displ_max", "fraction", 
            "Maximum marker displacement as a fraction of cell size in meters"
            ),
        itype_move_markers = ParameterInt(
            4, "itype_move_markers", "None", 
            "Displacement Options: 0 no motion 1 simp. advection; 4 4-th order Runge-Kutta"
            ),
        stype_move_markers = ParameterStr(
            "None", "stype_move_markers", "None", "Displacement options name"),
        iuse_local_adaptive_time_stepping = ParameterInt(
            0, "iuse_local_adaptive_time_stepping", "None", 
            "Activate local adaptive time stepping where displacement is calculated at staggered grid nodes. "
            * "The default value is 0 (off). If deactivated then adaptive time stepping using average grid spacing."
            ),

        # Markers - Subgrid diffusion parameters
        subgrid_diff_coef_stress = ParameterFloat(
            1.0, "subgrid_diff_coef_stress", "None", "Numerical subgrid stress diffusion coefficient"),
        subgrid_diff_coef_temp = ParameterFloat(
            1.0, "subgrid_diff_coef_temp", "None", "Numerical subgrid temperature diffusion coefficient"),

        )
end

function get_materials_parameters()::NamedTuple
    return (
        # *****************************
        # Materials Parameters
        # *****************************
        
        # Materials - Stress limits powerlaw parameters
        powerlaw_stress_min = ParameterFloat(
            1e4, "powerlaw_stress_min", "Pa", "Minimum stress limit for power law"),

        # Materials - Hydrothermal parameters
        iuse_hydrothermal = ParameterInt(
            0, "iuse_hydrothermal", "None", 
            "Hydrothermal circulation: 0 off; 1 on (iuse_topo must be 1)"
            ),
        iuse_melt_lens = ParameterInt(
            0, "iuse_melt_lens", "None", "Include shallow melt lens effect: 0 off; 1 on"),
        hydrothermal_smoothing_factor = ParameterFloat(
            0.75, "hydrothermal_smoothing_factor", "None", "Hydrothermal smoothing factor"),
        hydrothermal_nusselt_number = ParameterFloat(
            4.0, "hydrothermal_nusselt_number", "None", "Hydrothermal Nusselt number"),
        hydrothermal_max_temperature = ParameterFloat(
            600.0, "hydrothermal_max_temperature", "C", "Hydrothermal max temperature in Celsius"),
        hydrothermal_max_depth = ParameterFloat(
            6000.0, "hydrothermal_max_depth", "m", "Hydrothermal max depth in meters"),
        iuse_plastic_strain_rate_for_hydrothermal = ParameterInt(
            0, "iuse_plastic_strain_rate_for_hydrothermal", "None", 
            "Use plastic strain rate for hydrothermal circulation: 0 off; 1 on"
            ),
        hydrothermal_decay_length = ParameterFloat(
            50000.0, "hydrothermal_decay_length", "m", "Hydrothermal decay length in meters"),
        hydrothermal_buffer_distance = ParameterFloat(
            0.0, "hydrothermal_buffer_distance", "m", "Hydrothermal buffer distance in meters"),
        hydrothermal_plastic_strain_rate_reference = ParameterFloat(
            1e-13, "hydrothermal_plastic_strain_rate_reference", "1/s", 
            "Hydrothermal plastic strain rate reference in 1/s"
            ),
        iuse_plastic_strain_for_hydrothermal = ParameterInt(
            0, "iuse_plastic_strain_for_hydrothermal", "None", 
            "Use plastic strain for hydrothermal circulation: 0 off; 1 on"
            ),
        hydrothermal_plastic_strain_reference = ParameterFloat(
            0.5, "hydrothermal_plastic_strain_reference", "None", 
            "Hydrothermal plastic strain reference"
            ),
        sediment_thickness_threshold = ParameterFloat(
            2500.0, "sediment_thickness_threshold", "m", "Sediment thickness threshold in meters"),

        # Materials - Melt damage parameters
        iuse_melt_damage = ParameterInt(
            0, "iuse_melt_damage", "None", 
            "Integer flag that activates melt damage model where the damage "
            *"factor is calculated for markers in the melt damage zone above "
            *"partially molten mantle."
            ),
        melt_damage_distance = ParameterFloat(
            5000.0, "melt_damage_distance", "meters", 
            "Distance in meters from the shallow partially molten mantle where the damage "
            *"factor is calculated for markers in the melt damage zone."
            ),
        melt_damage_factor = ParameterFloat(
            10.0, "melt_damage_factor", "None", 
            "Maximum melt damage factor. Friction angle is divided by this factor "
            *"to account for melt damage."
            ),
        melt_damage_taper_distance = ParameterFloat(
            1000.0, "melt_damage_taper_distance", "meters", 
            "Distance in meters over which melt damage factor tapers to one."
            ),
        iuse_probabilistic_melt_damage = ParameterInt(
            0, "iuse_probabilistic_melt_damage", "None", 
            "Integer flag that activates probabilistic melt damage model where "
            *"the damage factor is calculated probabilistically for markers in "
            *"the melt damage zone."
            ),
        maximum_damage_probability = ParameterFloat(
            0.75, "maximum_damage_probability", "fraction", 
            "Maximum probability (fraction) of melt damage for markers in the central melt damage zone."
            ),
        intermediate_damage_probability = ParameterFloat(
            0.3, "intermediate_damage_probability", "fraction", 
            "Intermediate probability (fraction) of melt damage for markers in the central melt damage zone."
            ),
        magmatic_crust_height_threshold = ParameterFloat(
            1000.0, "magmatic_crust_height_threshold", "meters", 
            "Threshold height in meters for magmatic crust used in linear melt damage probability model"
            ),
        magmatic_crust_height_minimum = ParameterFloat(
            3000.0, "magmatic_crust_height_minimum", "meters", 
            "Minimum height in meters for magmatic crust used in linear melt damage probability model"
            ),
        magmatic_crust_height_maximum = ParameterFloat(
            10000.0, "magmatic_crust_height_maximum", "meters", 
            "Maximum height in meters for magmatic crust used in linear melt damage probability model."
            ),
        magmatic_crust_height_intermediate = ParameterFloat(
            2000.0, "magmatic_crust_height_intermediate", "meters", 
            "Intermediate height in meters for magmatic crust used in linear melt damage probability model."
            ),

        # Materials - Serpentinization parameters
        iuse_serpentinization = ParameterInt(
            0, "iuse_serpentinization", "None", 
            "Serpentinization model: 0 off; 1 on. Topography must be activated for serpentinization "
            *"(i.e. `iuse_topo` must be equal to 1). See [EarthBox.SurfaceProcesses.Topography.initialize!](@ref)) "
            *"for more details on topography initialization."
            ),
        serpentinization_temperature = ParameterFloat(
            600.0, "serpentinization_temperature", "C", 
            "Temperature in Celsius that controls when serpentinization rate reaches maximum and then rapidly "
            *"decreases. Typical values from the literature are around 613.15K (340.15C) "
            *"(e.g. [mezri24](@citet))."
            ),
        maximum_serpentinization_depth = ParameterFloat(
            20000.0, "maximum_serpentinization_depth", "m", 
            "Maximum submud serpentinization depth in meters"
            ),
        maximum_serpentinization_rate = ParameterFloat(
            1e-10, "maximum_serpentinization_rate", "1/s", 
            "Maximum serpentinization rate (1/s). Typical values from the literature are around "
            *"1e-10 to 1e-11 1/s (e.g. [mezri24](@citet))."
            ),
        nominal_strain_rate_serpentinization = ParameterFloat(
            1e-13, "nominal_strain_rate_serpentinization", "1/s", 
            "Nominal plastic strain rate for serpentinization (1/s). Plastic strain rate beyond "
            *"which the effect of plastic strain rate on serpentinization rate rapidly increases "
            *"and then becomes constant (e.g. Merzi et al., 2024). Typical values from the "
            *"literature are around 1e-12 to 1e-13 1/s (e.g. [mezri24](@citet))."
            ),
        # Materials - Description parameters
        itype_mat = ParameterInt(
            0, "itype_mat", "None", "Material and marker initialization model option"),
        stype_mat = ParameterStr(
            "None", "stype_mat", "None", "Material and marker initialization model option name"),
        itype_plasticity = ParameterInt(
            0, "itype_plasticity", "None", "Plasticity model type: 0 = viscoelastic; 1 = purely elastic"),
        stype_plasticity = ParameterStr(
            "None", "stype_plasticity", "None", "Plasticity model name"),

        # Materials - Random friction parameters
        iuse_random_fric = ParameterInt(
            0, "iuse_random_fric", "None", 
            "Randomize initial friction coefficient using delta_fric_coef: 0 off; 1 on"
            ),
        delta_fric_coef = ParameterFloat(
            0.1, "delta_fric_coef", "None", 
            "Friction coefficient perturbation used to randomize initial marker friction coefficients"
        ),
        iuse_random_fric_time = ParameterInt(
            0, "iuse_random_fric_time", "None", 
            "Randomize marker friction coefficients with each time step: 0 off; 1 on"
            ),
        randomization_factor = ParameterFloat(
            10.0, "randomization_factor", "None", 
            "Randomization factor used to randomize marker friction coefficients with each time step"
            ),

        # Materials - Stress limits yield parameters
        yield_stress_min = ParameterFloat(
            0.0, "yield_stress_min", "Pa", "Minimum stress limit for plastic failure in Pa"),
        yield_stress_max = ParameterFloat(
            1e38, "yield_stress_max", "Pa", "Maximum stress limit for plastic failure in Pa"),
        iuse_fluid_pressure_for_yield = ParameterInt(
            0, "iuse_fluid_pressure_for_yield", "None", 
            "Flag to use fluid pressure for yield stress calculation"
            ),
        plastic_healing_rate = ParameterFloat(
            0.0, "plastic_healing_rate", "1/s", "Rate of plastic healing in 1/s"),

        # Materials - Boundary friction parameters
        boundary_friction_width = ParameterFloat(
            0.0, "boundary_friction_width", "m", "Width of boundary friction zone"),
        boundary_friction_angle = ParameterFloat(
            0.0, "boundary_friction_angle", "Degrees", "Friction angle of boundary friction zone"),
        boundary_cohesion = ParameterFloat(
            0.0, "boundary_cohesion", "Pa", "Cohesion of boundary friction zone"),

        # Materials - Compaction parameters
        iuse_sed_porosity = ParameterInt(
            0, "iuse_sed_porosity", "None", "Depth-dependent sediment porosity: 0 off; 1 on"),
        conductivity_water = ParameterFloat(
            0.61, "conductivity_water", "W/m/K", "Thermal conductivity of water"),
        density_water = ParameterFloat(
            1000.0, "density_water", "kg/m^3", "Density of water"),
        heat_capacity_water = ParameterFloat(
            3000.0, "heat_capacity_water", "J/K/kg", "Heat capacity of water"),
        porosity_at_mudline = ParameterFloat(
            0.0, "porosity_at_mudline", "fraction", "Initial porosity at the sediment-water interface"),

        # Materials - Softening parameters
        iuse_viscous_strain_soft = ParameterInt(
            0, "iuse_viscous_strain_soft", "None", "Viscous strain softening: 0 off; 1 on"),
        vsoftfac = ParameterFloat(
            30.0, "vsoftfac", "None", "Softening factor multiplied by power-law pre-exponential"),

        # Materials - Viscosity limits parameters
        viscosity_min = ParameterFloat(
            1e18, "viscosity_min", "Pa.s", "Minimum viscosity"),
        viscosity_max = ParameterFloat(
            1e25, "viscosity_max", "Pa.s", "Maximum viscosity"),

        )
end

function get_melting_parameters()::NamedTuple
    return (
        # *****************************
        # Melting Parameters
        # *****************************
        
        # Melting - Rheology parameters
        viscosity_melt = ParameterFloat(
            1e17, "viscosity_melt", "Pa.s", "Viscosity of molten rock in Pa.s"),

        # Melting - Extrusion parameters
        iuse_extrusion = ParameterInt(
            0, "iuse_extrusion", "None", 
            "Integer flag that activates melt extrusion model: 0 off; 1 on"
            ),
        extrusion_volume_factor = ParameterFloat(
            0.2, "extrusion_volume_factor", "None",
            "Minimum value for the extrusion volume factor. "  
            *"The volume of lava produced at the surface of the model is calculated "
            *"as the product of the extrusion volume factor and the total volume "
            *"of magma produced during a time step. "
            *"The extrusion volume factor changes as a function of the characteristic "
            *"magmatic crust height. The characteristic magmatic crust height is "
            *"the height of new crust that would be formed if all "
            *"the melt in the mantle was extracted and emplaced at the top of the "
            *"model in a column with a width equal to the full extension velocity "
            *"times the time step. "
            *"This minimum value is used when the magmatic crust height is below the minimum "
            *"characteristic magmatic crust height. A linear model is used to define "
            *"the extrusion volume factor when the calculated magmatic crust height is between "
            *"the minimum and maximum characteristic magmatic crust height."
            ),
        extrusion_volume_factor_max = ParameterFloat(
            0.2, "extrusion_volume_factor_max", "None", 
            "Maximum value for the extrusion volume factor. This maximum "
            *"value is reached when the characteristic magmatic crust height is equal to the "
            *"maximum characteristic magmatic crust height."
            ),
        characteristic_magmatic_crust_height = ParameterFloat(
            0.0, "characteristic_magmatic_crust_height", "m", 
            "Characteristic height of magmatic crust in meters calculated based on "
            *"the total volume of magma produced during a given time step and the "
            *"full extensional velocity. The characteristic magmatic crust height "
            *"is the height of new crust that would be formed if all "
            *"the melt in the mantle was extracted and emplaced at the top of the "
            *"model in a column with a width equal to the full extension velocity "
            *"times the time step."
            ),
        characteristic_magmatic_crust_height_min = ParameterFloat(
            8000.0, "characteristic_magmatic_crust_height_min", "m", 
            "Minimum reference value for the characteristic magmatic crust height in meters. When "
            *"the calculated characteristic magmatic crust height is below this value the "
            *"extrusion volume factor is set to the minimum value "
            *"(i.e. `extrusion_volume_factor`)."
            ),
        characteristic_magmatic_crust_height_max = ParameterFloat(
            8000.0, "characteristic_magmatic_crust_height_max", "m", 
            "Maximum reference value for the characteristic magmatic crust height in meters. When "
            *"the calculated characteristic magmatic crust height is above this value the "
            *"extrusion volume factor is set to the maximum value "
            *"(i.e. `extrusion_volume_factor_max`)."
            ),
        characteristic_flow_length_subaerial = ParameterFloat(
            100000.0, "characteristic_flow_length_subaerial", "m", 
            "Characteristic flow length for subaerial eruptions in meters. Flow "
            *"distance of basaltic flows on land may be hundreds of kilometers (e.g. Ginkgo "
            *"flow of the Frenchman Springs Member of the Columbia River Basalts "
            *"has a flow length of 500 km)."
            ),
        characteristic_flow_length_submarine = ParameterFloat(
            10000.0, "characteristic_flow_length_submarine", "m", 
            "Characteristic flow length for submarine eruptions in meters."
            ),
        residual_lava_thickness_subaerial = ParameterFloat(
            30.0, "residual_lava_thickness_subaerial", "m", 
            "Residual lava thickness for subaerial eruptions in meters"
            ),
        residual_lava_thickness_submarine = ParameterFloat(
            60.0, "residual_lava_thickness_submarine", "m", 
            "Residual lava thickness for submarine eruptions in meters"
            ),
        iuse_random_eruption_location = ParameterInt(
            0, "iuse_random_eruption_location", "None", 
            "Integer flag that activates random eruption location: 0 off; 1 on"
            ),
        iuse_normal_eruption_location = ParameterInt(
            0, "iuse_normal_eruption_location", "None", 
            "Integer flag that activates normal distribution for eruption location: 0 off; 1 on"
            ),
        decimation_factor = ParameterInt(
            1, "decimation_factor", "None", 
            "Topography grid decimation factor for 1D lava flow model grid"
            ),
        initial_magma_flush_steps = ParameterInt(
            0, "initial_magma_flush_steps", "None", 
            "Number of time steps for initial magma flush"
            ),
        magma_flush_factor = ParameterFloat(
            0.8, "magma_flush_factor", "None", 
            "Factor that determines the amount of magma to flush for flush during magma flush time steps"
            ),
        width_eruption_domain_fixed = ParameterFloat(
            0.0, "width_eruption_domain_fixed", "m", 
            "Fixed width of the eruption domain in meters. If this is set to zero, the "
            *"width of the eruption domain is determined by the width of the "
            *"gabbroic molten zone including normal and layered gabbro."
            ),
        width_eruption_domain_fixed_max = ParameterFloat(
            0.0, "width_eruption_domain_fixed_max", "m", 
            "Maximum width of the eruption domain. If this is set to zero, the "
            *"width of the eruption domain is determined by the width of the "
            *"gabbroic molten zone including normal and layered gabbro."
            ),
        porosity_initial_lava_flow = ParameterFloat(
            0.0, "porosity_initial_lava_flow", "fraction", 
            "Initial porosity of lava flow in fraction used to decompact extruded gabbroic "
            *"magma to account for vesicles. This parameter is not used to "
            *"compact lava flows during burial. Instead, compaction properties "
            *"associated with the material type are used. This allows for the "
            *"decompaction of lava flows during eruption to account for porosity "
            *"from vesicles and the preservation of thickness during burial due to "
            *"chemical compaction and rigidity of the basalt matrix."
            ),
        decay_depth_lava_flow = ParameterFloat(
            2000.0, "decay_depth_lava_flow", "m", 
            "Depth at which the porosity of lava flow decays to zero in meters. This "
            *"parameter is used to compact lava flows during burial. Instead, "
            *"compaction properties associated with the material type are used. "
            *"This allows for the decompaction of lava flows during eruption to "
            *"account for porosity from vesicles and the preservation of thickness "
            *"during burial due to chemical compaction and rigidity of the basalt "
            *"matrix."
            ),
        time_of_next_eruption_myr = ParameterFloat(
            0.0, "time_of_next_eruption_myr", "Myr", "Time of next eruption in Myr"),
        eruption_interval_yr = ParameterFloat(
            100000.0, "eruption_interval_yr", "yr", "Eruption interval in years"),
        iuse_eruption_interval = ParameterInt(
            0, "iuse_eruption_interval", "None", 
            "Integer flag that activates eruption interval: 0 off; 1 on"
            ),

        # Melting - Extraction parameters
        melt_residual = ParameterFloat(
            0.0, "melt_residual", "None", "Sum of residual fractional melt particles"),
        ext_vol = ParameterFloat(
            0.0, "ext_vol", "m^3", "Volume of extruded material in m^3"),
        xmid_mol = ParameterFloat(
            0.0, "xmid_mol", "m", "x-location of the molten mantle zone in meters"),
        ytop_mol = ParameterFloat(
            0.0, "ytop_mol", "m", "y-location of the top of the molten mantle zone in meters"),
        width_mol = ParameterFloat(
            0.0, "width_mol", "m", "Width of the molten mantle zone in meters"),
        ndrainage_basin = ParameterInt(
            1, "ndrainage_basin", "None", "Number of drainage basins in the model"),
        smoothing_radius_drainage = ParameterFloat(
            2000.0, "smoothing_radius_drainage", "m", 
            "Smoothing radius (meters) used to smooth the top of the partially "
            *"molten mantle domain before calculating local maxima used to define "
            *"migration domains for melt extraction."
            ),
        characteristic_injection_width = ParameterFloat(
            2500.0, "characteristic_injection_width", "m", 
            "Characteristic width (meters) of magma injection into the shallow "
            *"mantle above the local maximum of the mantle partial melt domain. "
            *"This parameters is used only if use_shallow_mantle_injection is True. "
            *"Note that the actual injection width may be increased if the height "
            *"of injected magma exceeds the magma height limit."
            ),
        magma_height_limit = ParameterFloat(
            30000.0, "magma_height_limit", "m", 
            "Maximum height (meters) of magma that can be injected into the "
            *"mantle domain. This is used to limit the height of magma bodies "
            *"that are injected at the top of the mantle domain by increasing the "
            *"injection width."
            ),
        fractionation_threshold_limit = ParameterFloat(
            2000.0, "fractionation_threshold_limit", "m", 
            "Threshold limit (meters) for gabbroic fractionation. This is the "
            *"distance from Moho where the composition of gabbroic magma is changed "
            *"to fractionated or layered gabbroic magma to simulate the effects of "
            *"rapid fractionation in sills and the downward flow of the gabbro "
            *"glacier. This change in composition involves an increase in the solidus "
            *"leading to regions of partially molten fractioned gabbro as opposed "
            *"to large regions of pure gabbroic magma inconsistent with geophysical "
            *"observations."
            ),
        emplacement_temperature = ParameterFloat(
            1473.0, "emplacement_temperature", "K", 
            "Temperature at which magma is emplaced in Kelvin. "
            *"This parameter is currently not used in the code. Instead a default value of 1473.0 K is used."
            *"See function `transform_marker_to_magma` for details."
            ),
        number_of_injection_subdomains = ParameterInt(
            10, "number_of_injection_subdomains", "None", 
            "The number used to divide the characteristic injection width into "
            *"separate domains that are selected using a normal distribution if "
            *"`iuse_random_injection_subdomain = 1`."
            ),
        maximum_shallow_injection_depth = ParameterFloat(
            200000.0, "maximum_shallow_injection_depth", "m", 
            "Maximum submud depth (meters) for normal shallow mantle injection. "
            *"If the local maximum of the partially molten zone is deeper than "
            *"this value then the injection width is increased by a factor of 5. "
            *"THIS IS CURRENTLY TURNED OFF IN CODE. See function "
            *"extract_partial_melt_and_make_magma_body for details."
            ),
        extraction_fraction = ParameterFloat(
            1.0, "extraction_fraction", "fraction", 
            "Fraction of extracted melt that is extracted from the partially "
            *"molten zone. The remaining melt is left in the partially molten "
            *"zone and refertalizes the mantle."
            ),
        smoothing_radius_fractionation = ParameterFloat(
            10000.0, "smoothing_radius_fractionation", "m", 
            "Smoothing radius (meters) used to smooth the oceanic Moho for the  "
            *"gabbroic fractionation model whereby fractionation occurs based on "
            *"distance to Moho. Smoothing the Moho avoids spikes associated with "
            *"residual gabbroic particles in the mantle associated with magmatic "
            *"crust thinning."
            ),
        mantle_search_width = ParameterFloat(
            200000.0, "mantle_search_width", "m", 
            "Width (meters) of the mantle domain used to search for the shallowest "
            *"mantle marker in a subset of the model domain to improve computational "
            *"efficiency. This parameter is used only if use_shallow_mantle_injection "
            *"is True and may need to be adjusted depending on the characteristic "
            *"injection width. Adjustments should be made to ensure that the search "
            *"domain is equal to or larger than the injection domain."
            ),
        ndrainage_basin_old = ParameterInt(
            0, "ndrainage_basin_old", "None", "Old number of drainage basins"),

        # Melting - Options parameters
        iuse_melting = ParameterInt(
            0, "iuse_melting", "None", 
            "Integer flag that activates melt model with melt fraction calculation: 0 off; 1 on"
            ),
        iuse_melt_viscosity = ParameterInt(
            0, "iuse_melt_viscosity", "None", 
            "Integer flag that activates melt viscosity model where the viscosity of partially "
            *"molten rocks is set to melt viscosity if melt fraction exceeds 10%: 0 off; 1 on"
            ),
        iuse_melt_thermal_props = ParameterInt(
            0, "iuse_melt_thermal_props", "None", 
            "Integer flag that activates melt to impact thermal properties including density, heat "
            *"capacity in the term rho*cp and thermal expansivity in the adiabatic heating term alpha*T: 0 off; 1 on"
            ),
        iuse_extraction = ParameterInt(
            0, "iuse_extraction", "None", 
            "Integer flag that activates melt extraction model (0 = off, 1 = on). Extraction is required "
            *"for intrusive and extrusive processes."
            ),
        iuse_gabbroic_fractionation = ParameterInt(
            0, "iuse_gabbroic_fractionation", "None", 
            "Integer flag that activates gabbroic fractionation model: 0 off; 1 on. This "
            *"model is used to approximate the effects of rapid fractionation in "
            *"sills that produce layered gabbroic magma. This is accomplished by "
            *"transforming gabbroic magma to a fractionated gabbroic magma "
            *"composition once magma enters the crustal domain at the Moho."
            ),
        iuse_shallow_mantle_injection = ParameterInt(
            0, "iuse_shallow_mantle_injection", "None", 
            "Integer flag that activates shallow mantle injection model: 0 off; 1 on. This "
            *"model is used to inject mantle-derived magma at the shallowest part "
            *"of the mantle domain. If this option is False then mantle-derived "
            *"magma is injected at the shallowest location of the partially "
            *"molten mantle domain. This option approximates rapid ascent of melt "
            *"from the mantle to the crustal domain via dykes and channelized "
            *"melt networks."
            ),
        iuse_random_injection_subdomain = ParameterInt(
            0, "iuse_random_injection_subdomain", "None", 
            "Integer flag that activates the random injection model (0 = off, 1 = on) whereby "
            *"the characteristic injection width (characteristic_injection_width) is divided into "
            *"subdomains (number_of_injection_subdomains) that are selected using a uniform distribution. "
            *"The shallowest mantle marker within a selected subdomain is converted into magma. If this option "
            *"is False then the shallowest mantle marker within the injection width is converted into magma."
            ),
        iuse_normal_injection_subdomain = ParameterInt(
            0, "iuse_normal_injection_subdomain", "None", 
            "Integer flag that activates the normal injection model (0 = off, 1 = on) whereby the "
            *"characteristic injection width (characteristic_injection_width) is divided into "
            *"subdomains (number_of_injection_subdomains) that are selected using a normal distribution. "
            *"The shallowest mantle marker within a selected subdomain is converted into magma. If this option "
            *"is False then a uniform distribution is used to select the subdomain."
            ),
        iuse_depletion_density = ParameterInt(
            0, "iuse_depletion_density", "None", 
            "Integer flag that activate depletion density model (0 = off, 1 = on) whereby the "
            *"density of rocks that have undergone melt extraction is modified to account for "
            *"the partitioning of heavy elements into the melt."
            ), 
        iuse_exponential_viscosity_reduction = ParameterInt(
            0, "iuse_exponential_viscosity_reduction", "None", 
            "Integer flag that activates exponential viscosity reduction for partially "
            *"molten rocks whereby the effective viscosity of partially molten rocks is "
            *"reduced exponentially with extractable melt fraction for mantle rocks that are impacted by"
            *"melt extraction and melt fraction for crustal rocks that are currently not impacted "
            *"by melt extraction (0 = off, 1 = on). If this parameter is set to 0, the effective "
            *"viscosity of partially molten rocks is set to melt viscosity if melt fraction exceeds 10%."
            ),

        # Melting - Material IDs parameters
        ipmf1 = ParameterInt(
            5, "ipmf1", "None", "Material ID of the mantle"),
        ipmf2 = ParameterInt(
            6, "ipmf2", "None", "Material ID of the mantle"),
        ipmf3 = ParameterInt(
            7, "ipmf3", "None", "Material ID of the mantle"),
        ipmf4 = ParameterInt(
            8, "ipmf4", "None", "Material ID of the mantle"),
        ipmf_molten = ParameterInt(
            11, "ipmf_molten", "None", "Material ID of molten mantle"),
        ipmf_solid = ParameterInt(
            12, "ipmf_solid", "None", "Material ID of solidified mantle"),
        ipmf_molten_vol = ParameterInt(
            18, "ipmf_molten_vol", "None", "Material ID of molten extruded mantle"),
        ipmf_solid_vol = ParameterInt(
            19, "ipmf_solid_vol", "None", "Material ID of solidified extruded molten mantle"),

        )
end

function get_stokes_continuity_parameters()::NamedTuple
    return (
        # *****************************
        # Stokes-Continuity Parameters
        # *****************************
        
        # Stokes-continuity - Density parameters
        itype_density = ParameterInt(
            0, "itype_density", "None", "Integer ID of density model"),
        stype_density = ParameterStr(
            "Liao14", "stype_density", "None", "String name of density model"),

        # Stokes-continuity - Build parameters (requires constructor args - using defaults),
        ibuild_stokes = ParameterInt(
            1, "ibuild_stokes", "None", "0 define full system, 1 only define non-zero elements"),
        hshift_to_vxR = ParameterInt(
            300, "hshift_to_vxR", "None", 
            "Horizontal shift index for system building equal to (ynum-1)*3"
            ),
        Nstokes = ParameterInt(
            30000, "Nstokes", "None", 
            "Number of rows/columns in NxN matrix equal to (xnum-1)*(ynum-1)*3"
            ),
        nonzero_max_stokes = ParameterInt(
            316231, "nonzero_max_stokes", "None", 
            "Maximum number of non-zero values allowed set equal to xnum*ynum*31"
            ),
        pscale = ParameterFloat(
            1.0, "pscale", "None", "Coefficient for scaling pressure"),
        iuse_interface_stabilization = ParameterInt(
            0, "iuse_interface_stabilization", "None", 
            "Flag for interface stabilization: 0 off; 1 on"
            ),

        # Stokes-continuity - Solution norms parameters
        dsoluv1_abs_inf = ParameterFloat(
            1e38, "dsoluv1_abs_inf", "None", "Infinity norm of solution change"),
        dsoluv1_rel_inf = ParameterFloat(
            1e38, "dsoluv1_rel_inf", "None", "Relative infinity norm of solution change"),
        dsoluv1_abs_L2 = ParameterFloat(
            1e38, "dsoluv1_abs_L2", "None", "L2 norm of solution change"),
        dsoluv1_rel_L2 = ParameterFloat(
            1e38, "dsoluv1_rel_L2", "None", "Relative L2 norm of solution change"),
        dvx1_abs_L2 = ParameterFloat(
            1e38, "dvx1_abs_L2", "None", "L2 norm of vx change"),
        dvx1_rel_L2 = ParameterFloat(
            1e38, "dvx1_rel_L2", "None", "Relative L2 norm of vx change"),
        dvy1_abs_L2 = ParameterFloat(
            1e38, "dvy1_abs_L2", "None", "L2 norm of vy change"),
        dvy1_rel_L2 = ParameterFloat(
            1e38, "dvy1_rel_L2", "None", "Relative L2 norm of vy change"),
        dpr1_abs_L2 = ParameterFloat(
            1e38, "dpr1_abs_L2", "None", "L2 norm of pressure change"),
        dpr1_rel_L2 = ParameterFloat(
            1e38, "dpr1_rel_L2", "None", "Relative L2 norm of pressure change"),
        dvxy_abs_inf = ParameterFloat(
            1e38, "dvxy_abs_inf", "None", "Infinity norm of velocity solution"),
        dvxy_rel_inf = ParameterFloat(
            1e38, "dvxy_rel_inf", "None", "Relative infinity norm of velocity solution"),
        dvxy_abs_L2 = ParameterFloat(
            1e38, "dvxy_abs_L2", "None", "L2 norm of velocity solution"),
        dvxy_rel_L2 = ParameterFloat(
            1e38, "dvxy_rel_L2", "None", "Relative L2 norm of velocity solution"),
        global_yield_error = ParameterFloat(
            1e38, "global_yield_error", "None", "Global yield stress error"),

        # Stokes-continuity - Residual norms parameters
        resnl_L2_ini = ParameterFloat(
            1e38, "resnl_L2_ini", "None", "Initial L2 norm of Stokes system"),
        resnl_L2 = ParameterFloat(
            1e38, "resnl_L2", "None", "Current L2 norm of Stokes system"),
        resnl_rel_L2 = ParameterFloat(
            1e38, "resnl_rel_L2", "None", "Relative L2 norm of Stokes system"),
        resx_L2 = ParameterFloat(
            1e38, "resx_L2", "None", "L2 norm of vx residual"),
        resy_L2 = ParameterFloat(
            1e38, "resy_L2", "None", "L2 norm of vy residual"),
        resc_L2 = ParameterFloat(
            1e38, "resc_L2", "None", "L2 norm of pressure residual"),
        resx_L2_ini = ParameterFloat(
            1e38, "resx_L2_ini", "None", "Initial L2 norm of vx residual"),
        resy_L2_ini = ParameterFloat(
            1e38, "resy_L2_ini", "None", "Initial L2 norm of vy residual"),
        resc_L2_ini = ParameterFloat(
            1e38, "resc_L2_ini", "None", "Initial L2 norm of pressure residual"),
        resnlx_L2 = ParameterFloat(
            1e38, "resnlx_L2", "None", "Non-linear L2 norm of vx residual"),
        resnly_L2 = ParameterFloat(
            1e38, "resnly_L2", "None", "Non-linear L2 norm of vy residual"),
        resnlc_L2 = ParameterFloat(
            1e38, "resnlc_L2", "None", "Non-linear L2 norm of pressure residual"),
        resnlx_rel_L2 = ParameterFloat(
            1e38, "resnlx_rel_L2", "None", "Relative non-linear L2 norm of vx residual"),
        resnly_rel_L2 = ParameterFloat(
            1e38, "resnly_rel_L2", "None", "Relative non-linear L2 norm of vy residual"),
        resnlc_rel_L2 = ParameterFloat(
            1e38, "resnlc_rel_L2", "None", "Relative non-linear L2 norm of pressure residual"),
        resnlx_L2_ini = ParameterFloat(
            1e38, "resnlx_L2_ini", "None", "Initial non-linear L2 norm of vx residual"),
        resnly_L2_ini = ParameterFloat(
            1e38, "resnly_L2_ini", "None", "Initial non-linear L2 norm of vy residual"),
        resnlc_L2_ini = ParameterFloat(
            1e38, "resnlc_L2_ini", "None", "Initial non-linear L2 norm of pressure residual"),

        # Stokes-continuity - Picard parameters
        nglobal = ParameterInt(
            1, "nglobal", "None", 
            "Maximum number of global plasticity iterations (a.k.a. Picard iterations)"
            ),
        tolerance_picard = ParameterFloat(
            1e-2, "tolerance_picard", "None", 
            "Convergence criterion for the global plasticity iterations (a.k.a. Picard iterations)"
            ),
        iconverge = ParameterInt(
            0, "iconverge", "None", "Convergence flag"),
        iglobal = ParameterInt(
            0, "iglobal", "None", "Iteration counter for global stokes loop"),
        itype_global = ParameterInt(
            1, "itype_global", "None", "Type of global loop: 0 = marker-base; 1 = node-based"),
        stype_global = ParameterStr(
            "None", "stype_global", "None", "Global loop option name"),

        # Stokes-continuity - Velocity calc options parameters
        itype_velocity = ParameterInt(0, "itype_velocity", "None", 
            "Velocity calculation options: -1 all zeros; 0 solve equations; 1 solid body rotation"),
        stype_velocity = ParameterStr(
            "None", "stype_velocity", "None", "Velocity calculation option name"),
        # Stokes-continuity - Output parameters
        outtest = ParameterInt(
            0, "outtest", "None", "0 off; 1 output grid, vis., RHS, and bc as text files"),
        outtest2 = ParameterInt(
            0, "outtest2", "None", "1 output sparse matrix, rhs and solution vector as text files"),
        outtest3 = ParameterInt(
            0, "outtest3", "None", "1 output sparse matrix, rhs and solution vector as npy files"),

        )
end

function get_timestep_parameters()::NamedTuple
    return (
        # *****************************
        # Timestep Parameters
        # *****************************
        
        # Timestep - Boundary displacement stopping parameters
        iuse_boundary_displacement = ParameterInt(0, "iuse_boundary_displacement", "None", 
            "Flag to use boundary displacement stopping: 0 = off; 1 = on. When active, the model will "
            *"stop when the boundary displacement reaches the limit specified by `displ_limit`"),
        displ_limit = ParameterFloat(
            500.0, "displ_limit", "m", "Displacement limit for stopping the model in meters"),
        iuse_extensional_strain = ParameterInt(
            0, "iuse_extensional_strain", "None",
            "Flag to use extensional strain stopping: 0 = off; 1 = on. When active, the model will "
            *"stop when the extensional strain reaches the limit specified by `strain_limit`"
            ),
        strain_limit = ParameterFloat(
            1.0, "strain_limit", "None", "Strain limit for stopping the model"),

        # Timestep - Single increase parameters
        iuse_single_timestep_increase = ParameterInt(
            0, "iuse_single_timestep_increase", "None", 
            "Flag to use single timestep increase: 0 = off; 1 = on. When active, the time step will be "
            *"increased once during the simulation at the time step specified by ntime_increase"
            ),
        ntime_increase = ParameterInt(100, "ntime_increase", "None", 
            "Time step number at which to increase the timestep. This is only used if "
            *"`iuse_single_timestep_increase` is 1"
            ),

        # Timestep - Main time loop parameters
        model_duration_myr = ParameterFloat(
            0.0, "model_duration_myr", "Myr", 
            "Model duration in Myr used to estimate the maximum number of time steps for extensional "
            * "model with symmetric extension. "
            * "If adaptive time stepping is used during the simulation due to displacements that exceed "
            * "the user specified cell fraction then the actual model duration may be shorter than "
            * "`model_duration_myr`. "
            * "A buffer time of 3.0 Myr is added to the model duration when calculating the maximum "
            * "number of time steps assuming that there will be some time step reduction. "
            * "However, the user may need to increase the `model_duration_myr` to achieve the "
            * "desired model duration. If you want full control over the number of time steps then "
            * "set `ntimestep_max` and `timestep_viscoelastic` directly when using `run_time_steps` "
            * "function."
            ),
        ntimestep_max = ParameterInt(
            -9999, "ntimestep_max", "None", "Maximum number of model time steps"),
        timesum = ParameterFloat(
            0.0, "timesum", "s", "Total model time in seconds in seconds"),
        ntimestep = ParameterInt(
            0, "ntimestep", "None", "Model time step counter"),
        timestep_viscoelastic = ParameterFloat(
            NaN, "timestep_viscoelastic", "s", 
            "The viscoelastic time step in seconds used in the viscoelastic formulation of the "
            *"Stokes-continuity equations. The main model time step (`timestep`) is set equal to the "
            *"viscoelastic time step at the start of each time loop iteration prior to defining the "
            *"right-hand side terms that include the time-dependent viscoelastic factor. The "
            *"viscoelastic time step may be modified due to changes in boundary conditions depending "
            *"on user defined options"
            ),
        timestep_viscoelastic_step1 = ParameterFloat(
            NaN, "timestep_viscoelastic_step1", "s", 
            "The viscoelastic time step in seconds for the first velocity step in seconds"
            ),
        timestep_viscoelastic_step2 = ParameterFloat(
            NaN, "timestep_viscoelastic_step2", "s", 
            "The viscoelastic time step in seconds for the second velocity step in seconds"
            ),
        timestep = ParameterFloat(
            10000.0, "timestep", "s", 
            "Main model time step in seconds used for forecasting viscoelastic stress, solving heat "
            *"conduction equation, advecting markers and advancing model time. The main model time step "
            *"may be reduced after solving the Stokes-continuity equations to ensure that markers do not "
            *"advect beyond a user defined limit. The model time step is reset to the viscoelastic time "
            *"step at the start of each time loop iteration. The heat equation loop may use time steps "
            *"that are smaller than the main model time step if temperature changes are too large"
            ),
        iupdate_timestep = ParameterInt(
            0, "iupdate_timestep", "None", 
            "0 = off; 1 = update timestep using maximum displacement"),

        # Timestep - Multiple increase parameters
        iuse_multiple_timestep_increase = ParameterInt(
            0, "iuse_multiple_timestep_increase", "None", 
            "Flag to use multiple viscoelastic time step increases: 0 = off; 1 = on. When active the "
            *"viscoelastic time step will be increased multiple times during the simulation at the time "
            *"steps specified by `ntime_increase_1`, `ntime_increase_2`, etc"),
        ntime_increase_1 = ParameterInt(
            100, "ntime_increase_1", "None", "First time step number at which to increase the timestep"),
        ntime_increase_2 = ParameterInt(
            100, "ntime_increase_2", "None", "Second time step number at which to increase the timestep"),
        ntime_increase_3 = ParameterInt(
            100, "ntime_increase_3", "None", "Third time step number at which to increase the timestep"),
        ntime_increase_4 = ParameterInt(
            100, "ntime_increase_4", "None", "Fourth time step number at which to increase the timestep"),
        time_increase_factor = ParameterFloat(
            2.0, "time_increase_factor", "None", "Multiplication factor used to increase the "
            *"viscoelastic time step"
            ),
        cell_displ_factor = ParameterFloat(
            2.0, "cell_displ_factor", "None", 
            "Factor used to adjust the maximum allowed cell displacement fraction during multiple "
            *"viscoelastic time step increases"
            ),

        # Timestep - Thermal loop parameters
        timestep_heat = ParameterFloat(
            0.0, "timestep_heat", "s", "Thermal time step in seconds used in adaptive heat solver loop"),
        timestep_sum = ParameterFloat(
            0.0, "timestep_sum", "s", "Total thermal loop time in seconds"),

        # Timestep - Output steps parameters
        timestep_out = ParameterFloat(
            1.0e3, "timestep_out", "s", "Time step for output"),
        nskip = ParameterInt(
            0, "nskip", "None", 
            "Number of time steps to skip before output is generated. This is calculated as "
            *"`floor(Int, timestep_out/timestep_viscoelastic)`"
            ),
        icount_output = ParameterInt(
            0, "icount_output", "None", 
            "Output loop counter used to track time steps between outputs. An output event is triggered "
            *"when this counter is equal to nskip, which involves setting the output counter to zero. "
            *"This counter is incremented by one at the end of each time step"
            ),
        noutput = ParameterInt(
            0, "noutput", "None", "Number of outputs generated during a model run"),

        )
end

function get_topography_parameters()::NamedTuple
    return (
        # *****************************
        # Topography Parameters
        # *****************************
        
        # Topography - Reference lithosphere parameters
        thickness_upper_continental_crust_ref = ParameterFloat(
            22000.0, "thickness_upper_continental_crust_ref", "m", 
            "Thickness of upper continental crust of the reference lithosphere in meters"
            ),
        thickness_lower_continental_crust_ref = ParameterFloat(
            10000.0, "thickness_lower_continental_crust_ref", "m", 
            "Thickness of lower continental crust of the reference lithosphere in meters"
            ),
        thickness_lithosphere_ref = ParameterFloat(
            125_000.0, "thickness_lithosphere_ref", "m", 
            "Thickness of lithosphere for the reference lithosphere in meters"
            ),
        gridy_spacing_ref = ParameterFloat(
            100.0, "gridy_spacing_ref", "m", 
            "Spacing of grid cells in y-direction for the reference lithosphere in meters"
            ),
        temperature_top_ref = ParameterFloat(
            0.0, "temperature_top_ref", "C", 
            "Temperature at top of lithosphere for the reference lithosphere in Celsius"
            ),
        temperature_moho_ref = ParameterFloat(
            600.0, "temperature_moho_ref", "C", 
            "Temperature at Moho of the reference lithosphere in Celsius"
            ),
        temperature_base_lith_ref = ParameterFloat(
            1330.0, "temperature_base_lith_ref", "C", 
            "Temperature at base of lithosphere for the reference lithosphere in Celsius"
            ),
        adiabatic_gradient_ref = ParameterFloat(
            0.4, "adiabatic_gradient_ref", "K/km", 
            "Adiabatic gradient below the reference lithosphere in K/km"
            ),
        iuse_linear_segments = ParameterInt(
            1, "iuse_linear_segments", "None", 
            "Integer flag that controls the type of temperature profile used in "
            *"the reference lithosphere model. If `iuse_linear_segments = 1`, a temperature profile "
            *"with four linear segments is used that is controlled by the following user defined "
            *"parameters: `temperature_top_ref`, `temperature_moho_ref`, and `temperature_base_lith_ref`. "
            *"If `iuse_linear_segments = 0`, the temperature profile is calculated using an analytical 3-layer model "
            *"defined during marker temperature. See the `AnalyticalThreeLayer` model description at "
            *"[MarkerTemperature](@ref)). Default is 1."
            ),

        # Topography - Downhill diffusion parameters
        downhill_diff_elev_max = ParameterFloat(
            1.0, "downhill_diff_elev_max", "m", "d2Yt - max elevation"),
        transport_length = ParameterFloat(
            1.0, "transport_length", "m", "Transport length scale in meters"),
        topo_diff_coef = ParameterFloat(
            0.0, "topo_diff_coef", "m^2/s", "Topography diffusion coefficient in m^2/s"),
        subaerial_slope_diffusivity = ParameterFloat(
            0.0, "subaerial_slope_diffusivity", "m^2/s", 
            "Subaerial slope diffusivity in m^2/s. Typical value used in the literature is 0.25 "
            *"m^2/yr (7.9e-9 m^2/s) (e.g. Andres-Martinez et al., 2019; Armitage et al., 2015)"
            ),
        precipitation_rate = ParameterFloat(0.0, "precipitation_rate", "m/s", 
            "Precipitation rate in m/s used to calculate water flux in drainage basins for "
            *"fluvial transport diffusivity model. Precipitation rate (m/s). Used to "
            *"calculate water flux in drainage basins. Typical value used in the literature is "
            *"1 m/yr (3.2e-8 m/s) (e.g. Andres-Martinez et al.; 2019; Huffman et al., 2009)."
            ),
        subaerial_transport_coefficient = ParameterFloat(
            0.0, "subaerial_transport_coefficient", "None", 
            "Subaerial discharge transport coefficient. Used to calculate effective subaerial fluvial "
            *"transport diffusivity that includes slope diffusivity, precipitation rate and downstream "
            *"distances. Typical values used in the literature are 1e-4 (low transport) to 1e-2 (high "
            *"transport) (e.g. Andres-Martinez et al., 2019; Armitage et al., 2015)."
            ),
        submarine_slope_diffusivity = ParameterFloat(
            0.0, "submarine_slope_diffusivity", "m^2/s", 
            "Maximum submarine slope diffusivity in m^2/s used in diffusivity model that decays exponentially "
            *"with water depth. Typical value used in the literature is 100 m^2/yr (3.2e-9 m^2/s) (e.g. "
            *"Andres-Martinez et al., 2019; Kaufman et al., 1991)."
            ),
        submarine_diffusion_decay_depth = ParameterFloat(
            0.0, "submarine_diffusion_decay_depth", "m", 
            "Submarine diffusion decay depth (m). Typical value used in the "
            *"literature is 1000-2000 m (e.g. Andres-Martinez et al., 2019; Kaufman"
            *"et al., 1991; Perez-Gussinye et al., 2020)."
            ),
        number_of_transport_timesteps_per_model_timestep = ParameterInt(
            5, "number_of_transport_timesteps_per_model_timestep", "None", 
            "Number of transport time steps per model time step."
            ),
        transport_timestep = ParameterFloat(
            NaN, "transport_timestep", "s", 
            "Transport timestep (s). Typical value used in the literature is "
            *"1000 years (3.1536e10 s) (e.g. Andres-Martinez et al., 2019)."
            ),
        iuse_compaction_correction = ParameterInt(
            0, "iuse_compaction_correction", "None", 
            "Integer flag activating compaction correction during sediment transport where newly "
            *"deposited sediment thickness and pre-existing sediment thickness "
            *"are corrected for compaction."
            ),

        # Topography - Depo and erosion rates parameters
        erosion_rate = ParameterFloat(
            0.0, "erosion_rate", "m/s", 
            "Erosion rate above water level in m/s"
            ),
        sedimentation_rate = ParameterFloat(
            0.0, "sedimentation_rate", "m/s", 
            "Sedimentation rate below water level in m/s"
            ),
        pelagic_sedimentation_rate = ParameterFloat(
            0.0, "pelagic_sedimentation_rate", "m/s", 
            "Pelagic sedimentation rate (m/s). Typical value used in the "
            *"literature are 0.3 mm/yr (syn-rift) to 0.01 mm/yr (post-rift) "
            *"(e.g. Perez-Gussinye et al., 2020)"
            ),
        pelagic_sedimentation_rate_reduction_factor = ParameterFloat(
            1.0, "pelagic_sedimentation_rate_reduction_factor", "None", 
            "The pelagic sedimentation rate is divided by this factor after "
            *"the specified pelagic sedimentation rate reduction time."
            ),
        pelagic_sedimentation_rate_reduction_time = ParameterFloat(
            5000.0, "pelagic_sedimentation_rate_reduction_time", "Myr", 
            "Time in Myr after which the pelagic sedimentation rate is reduced "
            *"by the pelagic_sedimentation_rate_reduction_factor."
            ),
        salt_deposition_rate = ParameterFloat(
            0.0, "salt_deposition_rate", "m/s", "Salt deposition rate in m/s"),
        salt_start_time = ParameterFloat(
            1000.0, "salt_start_time", "Myr", "Salt start time in Myr"),
        salt_end_time = ParameterFloat(
            1001.0, "salt_end_time", "Myr", "Salt end time in Myr"),

        # Topography - Sea level parameters
        itype_sealevel = ParameterInt(
            0, "itype_sealevel", "None", "Integer flag controlling the sealevel model"),
        stype_sealevel = ParameterStr(
            "None", "stype_sealevel", "None", "Sealevel option name"),
        y_water_ini = ParameterFloat(
            10000.0, "y_water_ini", "m", "Initial and current seawater level in meters"),
        base_level_shift = ParameterFloat(
            0.0, "base_level_shift", "m", 
            "Shift in base level in meters. Positive values shift the base level down in the +y direction."
            ),
        base_level_shift_end_time = ParameterFloat(
            -1.0, "base_level_shift_end_time", "Myr", "End time for base level shift in Myr"),
        y_sealevel = ParameterFloat(
            0.0, "y_sealevel", "m", "Seawater level y-coordinate in meters"),

        # Topography - Model options parameters
        iuse_topo = ParameterInt(
            0, "iuse_topo", "None", 
            "Integer flag used to enable topography model: 0 = off; 1 = on"
            ),
        itype_topo = ParameterInt(
            0, "itype_topo", "None", 
            "ID for topography node advection option: 0 = Runge Kutta with interpolation; 1 = antidiffusion"
            ),
        stype_topo = ParameterStr(
            "None", "stype_topo", "None", 
            "Topography advection option name"
            ),
        iuse_downhill_diffusion = ParameterInt(
            0, "iuse_downhill_diffusion", "None", 
            "Integer flag that activates sediment transport model with downhill diffusion "
            *"(0 = off, 1 = on). If option is set to off uniform erosion and deposition rates will be used"
            ),
        iuse_salt_deposition = ParameterInt(
            0, "iuse_salt_deposition", "None",
            "Integer flag that activates salt deposition model: 0 = off; 1 = on"
            ),

        # Topography - Topo grid parameters
        topo_xsize = ParameterFloat(
            100000.0, "topo_xsize", "m", "Topography model size in horizontal direction in meters"),
        toponum = ParameterInt(
            1001, "toponum", "None", "Number of nodes in the x-direction of the topography model"),
        dx_topo = ParameterFloat(
            NaN, "dx_topo", "m", "Topography grid step size in meters"),
        nsmooth_top_bottom = ParameterInt(
            2, "nsmooth_top_bottom", "None", 
            "Number of topography nodes that will be used to calculate a running average for the "
            *"calculated y-coordinate of the top and bottom of marker layers"
            ),
        marker_search_factor = ParameterFloat(
            2.0, "marker_search_factor", "None", 
            "Search factor used to calculate the search radius for determining the shallowest and "
            *"deepest markers associated with a given topography node"
            ),

        )
end

function get_material_model_parameters()::NamedTuple
    return (
        #********************
        # Material Parameters
        # *******************
        # Material parameters are defined in material input and library files.

        # Material Model Parameters
        name = ParameterStr("dummy", "name", "None", "Name of material"),
        ncolors = ParameterInt(0, "ncolors", "None", "Number of colors"),
        nmats = ParameterInt(
            100, "nmats", "None", "Maximum number of materials allowed in model"),
        material_model_id = ParameterInt(
            0, "material_model_id", "None", "Material model ID"),
        matid = ParameterInt(0, "matid", "None", "Material ID"),
        mat_name = ParameterStr("dummy", "mat_name", "None", "Material name"),
        mat_type = ParameterStr("dummy", "mat_type", "None", "Material type"),
        mat_domain = ParameterStr("dummy", "mat_domain", "None", "Material domain"),
        description = ParameterStr("dummy", "description", "None", "Material description"),

        # Density Parameters
        standard_density = ParameterFloat(
            0.0, "standard_density", "kg/m^3", "Standard density in kg/m^3"),
        thermal_expansion = ParameterFloat(
            0.0, "thermal_expansion", "1/K", "Thermal expansion coefficient in 1/K"),
        compressibility = ParameterFloat(
            0.0, "compressibility", "1/Pa", "Compressibility in 1/Pa"),
        melt_density = ParameterFloat(0.0, "melt_density", "kg/m^3", "Melt density in kg/m^3"),
        viscosity_iso = ParameterFloat(0.0, "viscosity_iso", "Pa.s", 
            "Isoviscous viscosity in Pa.s"),

        # Flow Law Parameters
        flow_type = ParameterInt(0, "flow_type", "None", "Flow law type"),
        flow_stype = ParameterStr("None", "flow_stype", "None", "Flow law type name"),

        # Dislocation Creep Parameters
        pre_exponential_dc = ParameterFloat(
            0.0, "pre_exponential_dc", "1/s/MPa^n", "Pre-exponential factor for dislocation creep in 1/s/MPa^n"),
        stress_exponent_n_dc = ParameterFloat(
            0.0, "stress_exponent_n_dc", "None", "Stress exponent for dislocation creep"),
        activation_energy_dc = ParameterFloat(
            0.0, "activation_energy_dc", "kJ/mol", "Activation energy for dislocation creep in kJ/mol"),
        activation_volume_dc = ParameterFloat(
            0.0, "activation_volume_dc", "J/MPa/mol", "Activation volume for dislocation creep in J/MPa/mol"),

        # Peierls Creep Parameters
        pre_exponential_pei = ParameterFloat(
            0.0, "pre_exponential_pei", "s^-m1*MPa^-m2", "Pre-exponential factor for Peierls creep in s^-m1*MPa^-m2"),
        stress_exponent_m1_pei = ParameterFloat(
            0.0, "stress_exponent_m1_pei", "None", "Stress exponent m1 for Peierls creep"),
        stress_exponent_m2_pei = ParameterFloat(
            0.0, "stress_exponent_m2_pei", "None", "Stress exponent m2 for Peierls creep"),
        peierls_stress = ParameterFloat(
            0.0, "peierls_stress", "MPa", "Peierls stress in MPa"),

        # Diffusion Creep Parameters
        pre_exponential_difc = ParameterFloat(
            0.0, "pre_exponential_difc", "1/s/MPa^n", "Pre-exponential factor for diffusion creep in 1/s/MPa^n"),
        activation_energy_difc = ParameterFloat(
            0.0, "activation_energy_difc", "kJ/mol", "Activation energy for diffusion creep in kJ/mol"),
        activation_volume_difc = ParameterFloat(
            0.0, "activation_volume_difc", "J/MPa/mol", "Activation volume for diffusion creep in J/MPa/mol"),

        # Temperature Dependent Viscosity Parameters
        pre_exponential_td = ParameterFloat(
            0.0, "pre_exponential_td", "Pa.s", 
            "Pre-exponential factor for temperature dependent viscosity in Pa.s"
            ),
        activation_energy_td = ParameterFloat(
            0.0, "activation_energy_td", "kJ/mol", 
            "Activation energy for temperature dependent viscosity in kJ/mol"
            ),

        # Blankenbach89 Viscosity Parameters
        viscosity_ref_blankenbach89 = ParameterFloat(
            0.0, "viscosity_ref_blankenbach89", "Pa.s", "Reference viscosity for Blankenbach89 in Pa.s"),
        b_term_blankenbach89 = ParameterFloat(
            0.0, "b_term_blankenbach89", "None", "b term for Blankenbach89"),
        c_term_blankenbach89 = ParameterFloat(
            0.0, "c_term_blankenbach89", "None", "c term for Blankenbach89"),

        # Shear Modulus Parameters
        shear_modulus = ParameterFloat(
            0.0, "shear_modulus", "Pa", "Shear modulus in Pa"),

        # Plasticity Parameters
        cohesion_initial = ParameterFloat(
            0.0, "cohesion_initial", "Pa", "Initial cohesion in Pa"),
        cohesion_final = ParameterFloat(
            0.0, "cohesion_final", "Pa", "Final cohesion in Pa"),
        friction_angle_initial = ParameterFloat(
            0.0, "friction_angle_initial", "Degrees", "Initial friction angle in degrees"),
        friction_angle_final = ParameterFloat(
            0.0, "friction_angle_final", "Degrees", "Final friction angle in degrees"),
        sine_friction_angle_initial = ParameterFloat(
            0.0, "sine_friction_angle_initial", "None", "Initial sine friction angle"),
        sine_friction_angle_final = ParameterFloat(
            0.0, "sine_friction_angle_final", "None", "Final sine friction angle"),
        strain_initial = ParameterFloat(
            0.0, "strain_initial", "None", "Initial strain"),
        strain_final = ParameterFloat(
            0.0, "strain_final", "None", "Final strain"),

        # Heat Capacity Parameters
        heat_capacity = ParameterFloat(
            0.0, "heat_capacity", "J/kg/K", "Heat capacity in J/kg/K"),

        # Thermal Conductivity Parameters
        thermal_conductivity_ref = ParameterFloat(
            0.0, "thermal_conductivity_ref", "W/m/K", "Reference thermal conductivity in W/m/K"),
        thermal_conductivity_a = ParameterFloat(
            0.0, "thermal_conductivity_a", "W/m", "Thermal conductivity parameter a in W/m"),

        # Radiogenic Heat Production Parameters
        radiogenic_heat_production = ParameterFloat(
            0.0, "radiogenic_heat_production", "W/m^3", "Radiogenic heat production in W/m^3"),

        # Material Color Parameters
        red_fraction = ParameterFloat(
            0.0, "red_fraction", "None", "Red color fraction"),
        green_fraction = ParameterFloat(
            0.0, "green_fraction", "None", "Green color fraction"),
        blue_fraction = ParameterFloat(
            0.0, "blue_fraction", "None", "Blue color fraction"),
        
        # Dilatation Angle Parameters
        dilatation_angle = ParameterFloat(
            0.0, "dilatation_angle", "Degrees", "Dilatation angle in degrees"),

        # Melting Parameters
        itype_solidus = ParameterInt(
            -1, "itype_solidus", "None", "Solidus type ID"),
        stype_solidus = ParameterStr(
            "None", "stype_solidus", "None", "Solidus type name"),
        itype_liquidus = ParameterInt(
            -1, "itype_liquidus", "None", "Liquidus type ID"),
        stype_liquidus = ParameterStr(
            "None", "stype_liquidus", "None", "Liquidus type name"),
        latent_heat = ParameterFloat(
            400000.0, "latent_heat", "J/kg", "Latent heat in J/kg"),

        # Compaction Parameters
        porosity_initial = ParameterFloat(
            0.0, "porosity_initial", "fraction", "Initial porosity in fraction"),
        porosity_decay_depth = ParameterFloat(
            2001.0, "porosity_decay_depth", "m", "Porosity decay depth in m"),

        )
end

function get_material_override_parameters()::NamedTuple
    return (
        #******************
        # Material Override
        # *****************
        # Material override parameters are used to override material properties that are defined in the 
        # material library file.

        # Material Override: Friction Angles
        friction_angle_initial_strong_zone = ParameterFloat(
            0.0, "friction_angle_initial_strong_zone", "degrees", 
            "Initial friction angle for strong zones"
            ),
        friction_angle_final_strong_zone = ParameterFloat(
            0.0, "friction_angle_final_strong_zone", "degrees", 
            "Final friction angle for strong zones"
            ),
        friction_angle_initial_solidified_basalt = ParameterFloat(
            0.0, "friction_angle_initial_solidified_basalt", "degrees", 
            "Initial friction angle for solidified basalt"
            ),
        friction_angle_final_solidified_basalt = ParameterFloat(
            0.0, "friction_angle_final_solidified_basalt", "degrees", 
            "Final friction angle for solidified basalt"
            ),
        friction_angle_initial_oceanic_crust = ParameterFloat(
            0.0, "friction_angle_initial_oceanic_crust", "degrees", 
            "Initial friction angle for oceanic crust"
            ),
        friction_angle_final_oceanic_crust = ParameterFloat(
            0.0, "friction_angle_final_oceanic_crust", "degrees", 
            "Final friction angle for oceanic crust"
            ),
        friction_angle_initial_sediment = ParameterFloat(
            0.0, "friction_angle_initial_sediment", "degrees", 
            "Initial friction angle for sediment"
            ),
        friction_angle_final_sediment = ParameterFloat(
            0.0, "friction_angle_final_sediment", "degrees", 
            "Final friction angle for sediment"
            ),

        # Material Override: Latent Heat and Solidus and Liquidus
        latent_heat_mantle = ParameterFloat(
            400000.0, "latent_heat_mantle", "J/kg", 
            "Latent heat for mantle"
            ),
        latent_heat_oceanic_crust = ParameterFloat(
            400000.0, "latent_heat_oceanic_crust", "J/kg", 
            "Latent heat for oceanic crust"
            ),
        mantle_solidus = ParameterStr(
            "PeridotiteKatz2003", "mantle_solidus", "None", 
            "String name of solidus model for mantle"
            ),
        mantle_liquidus = ParameterStr(
            "PeridotiteKatz2003", "mantle_liquidus", "None", 
            "String name of liquidus model for mantle"
            ),
        layered_gabbro_solidus = ParameterStr(
            "GabbroGlacier", "layered_gabbro_solidus", "None", 
            "String name of solidus model for layered gabbro"
            ),
        layered_gabbro_liquidus = ParameterStr(
            "GabbroGlacier", "layered_gabbro_liquidus", "None", 
            "String name of liquidus model for layered gabbro"
            ),
        gabbro_solidus = ParameterStr(
            "GabbroGerya2010", "gabbro_solidus", "None", 
            "String name of solidus model for gabbro"
            ),
        gabbro_liquidus = ParameterStr(
            "GabbroGerya2010", "gabbro_liquidus", "None", 
            "String name of liquidus model for gabbro"
            ),

        # Material Override: Pre-exponential factors
        scale_factor_crustal_dislocation_creep = ParameterFloat(
            1.0, "scale_factor_crustal_dislocation_creep", "None", 
            "Scale factor for pre-exponential term for crustal dislocation creep"
            ),
        scale_factor_mantle_dislocation_creep = ParameterFloat(
            1.0, "scale_factor_mantle_dislocation_creep", "None", 
            "Scale factor for pre-exponential term for mantle dislocation creep"
            ),
        scale_factor_mantle_diffusion_creep = ParameterFloat(
            1.0, "scale_factor_mantle_diffusion_creep", "None", 
            "Scale factor for pre-exponential term for mantle diffusion creep"
            ),
        scale_factor_oceanic_crust_dislocation_creep = ParameterFloat(
            1.0, "scale_factor_oceanic_crust_dislocation_creep", "None", 
            "Scale factor for pre-exponential term for oceanic crust dislocation creep"
            ),
        scale_factor_solidified_basalt_dislocation_creep = ParameterFloat(
            1.0, "scale_factor_solidified_basalt_dislocation_creep", "None", 
            "Scale factor for pre-exponential term for solidified basalt dislocation creep"
            ),
        scale_factor_sediment_dislocation_creep = ParameterFloat(
            1.0, "scale_factor_sediment_dislocation_creep", "None", 
            "Scale factor for pre-exponential term for sediment dislocation creep"
            ),

        # Material Override: Strain limits for weakening
        strain_initial_strong_zone = ParameterFloat(
            0.0, "strain_initial_strong_zone", "None", 
            "Initial reference strain used for linear strain weakening model for strong zones"
            ),
        strain_final_strong_zone = ParameterFloat(
            0.0, "strain_final_strong_zone", "None", 
            "Final reference strain used for linear strain weakening model for strong zones"
            ),
        strain_initial_basalt = ParameterFloat(
            0.0, "strain_initial_basalt", "None", 
            "Initial reference strain used for linear strain weakening model for solidified basalt"
            ),
        strain_final_solidified_basalt = ParameterFloat(
            0.0, "strain_final_solidified_basalt", "None", 
            "Final reference strain used for linear strain weakening model for solidified basalt"
            ),
        strain_initial_oceanic_crust = ParameterFloat(
            0.0, "strain_initial_oceanic_crust", "None", 
            "Initial reference strain used for linear strain weakening model for oceanic crust"
            ),
        strain_final_oceanic_crust = ParameterFloat(
            0.0, "strain_final_oceanic_crust", "None", 
            "Final reference strain used for linear strain weakening model for oceanic crust"
            ),
        strain_initial_sediment = ParameterFloat(
            0.0, "strain_initial_sediment", "None", 
            "Initial reference strain used for linear strain weakening model for sediment"
            ),
        strain_final_sediment = ParameterFloat(
            0.0, "strain_final_sediment", "None", 
            "Final reference strain used for linear strain weakening model for sediment"
            ),

        )
end

function get_other_parameters()::NamedTuple
    return (
        #*****************
        # Other Parameters
        #*****************
        iuse_mumps = ParameterInt(
            0, "iuse_mumps", "None", "Flag to use MUMPS solver: 0 off; 1 on"),
    )
end

end # module