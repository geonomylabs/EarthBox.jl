function get_heat_equation_parameters()::NamedTuple
    return @params (
        # *****************************
        # Heat Equation Parameters
        # *****************************
        
        # Heat equation - Temp change limit parameters
        max_temp_change = ParameterFloat(
            20.0, "K", "Maximal temperature change allowed for one timestep in Kelvin"),

        # Heat equation - Steady state parameters
        thick_thermal_lithosphere = ParameterFloat(
            NaN, "m", "Thickness of the thermal lithosphere in meters"),
        conductivity_upper_crust = ParameterFloat(
            3.0, "W/m/K", "Thermal conductivity of upper crust W/m/K"),
        conductivity_lower_crust = ParameterFloat(
            3.0, "W/m/K", "Thermal conductivity of lower crust W/m/K"),
        conductivity_mantle = ParameterFloat(
            3.0, "W/m/K", "Thermal conductivity of mantle W/m/K"),
        heat_production_upper_crust = ParameterFloat(
            1.5e-6, "W/m^3", "Radiogenic heat production in upper crust W/m^3"),
        heat_production_lower_crust = ParameterFloat(
            2.0e-7, "W/m^3", "Radiogenic heat production in lower crust W/m^3"),
        heat_production_mantle = ParameterFloat(
            0.0, "W/m^3", "Radiogenic heat production in mantle W/m^3"),
        amplitude_perturbation = ParameterFloat(
            0.0, "m", "Amplitude of central thermal perturbation in meters"),
        width_perturbation = ParameterFloat(
            0.0, "m", "Width of central thermal perturbation in meters"),
        temperature_surface = ParameterFloat(
            0.0, "K", 
            "Initial temperature at the surface for lithospheric models with sticky air in Kelvin"
            ),
        temperature_moho = ParameterFloat(
            0.0, "K", "Initial temperature at the Moho for lithospheric models in Kelvin"),
        temperature_base_lith = ParameterFloat(
            NaN, "K", "Initial temperature at the base of the lithosphere in Kelvin"),
        # This may be obsolete
        temperature_base_lithosphere = ParameterFloat(
            NaN, "K", 
            "Temperature at the base of the thermal lithosphere in Kelvin"
            ),

        # Heat equation - Initial condition parameters
        itype_temp = ParameterInt(
            0, "None", "Option flag for initial temperature condition"),
        stype_temp = ParameterStr(
            "None", "None", "Initial temperature condition option name"),
        temperature_uniform = ParameterFloat(
            0.0, "K", "Uniform initial temperature value in Kelvin"),
        temperature_of_wave = ParameterFloat(
            0.0, "K", "Temperature of wave in Kelvin"),
        age_lithosphere_left = ParameterFloat(
            0.0, "Myr", "Age of lithosphere to the left of fracture zone in Myr"),
        age_lithosphere_right = ParameterFloat(
            0.0, "Myr", "Age of lithosphere to the right of fracture zone in Myr"),
        thermal_lithosphere_depth_left = ParameterFloat(
            100000.0, "m", 
            "Thermal lithosphere depth to the left of fracture zone in meters"
            ),
        thermal_lithosphere_depth_right = ParameterFloat(
            200000.0, "m", 
            "Thermal lithosphere depth to the right of fracture zone in meters"
            ),
        thermal_diffusivity = ParameterFloat(
            1e-6, "m^2/s", "Thermal diffusivity in m^2/s"),
        adiabatic_gradient = ParameterFloat(
            0.4, "K/km", "Adiabatic gradient in K/km"),
        
        # Heat equation - RhoCp parameters
        itype_rhocp = ParameterInt(
            0, "None", "Integer flag for density*heat_capacity model"),
        stype_rhocp = ParameterStr(
            "None", "None", "Density*heat_capacity model option name"),
        maximum_heat_capacity = ParameterFloat(
            1200.0, "J/kg/K", 
            "Maximum heat capacity for rock properties used with Waples temperature-dependent model in J/kg/k"
            ),

        # Heat equation - Build parameters (requires constructor args - using defaults),
        ibuild_heat = ParameterInt(
            1, "None", "0 define full system, 1 only define non-zero elements"),
        Nheat = ParameterInt(
            10201, "None", "Number of rows or columns in the large NxN matrix"),
        nonzero_max_heat = ParameterInt(
            316231, "None", 
            "Maximum number of non-zero values allowed for the large matrix"
            ),

        # Heat equation - Heat options parameters
        iuse_heat = ParameterInt(
            1, "None", "Solve heat equation: 0 off; 1 on"),
        iuse_shear_heating = ParameterInt(
            0, "None", "Shear heating: 0 off; 1 on"),
        iuse_adiabatic_heating = ParameterInt(
            0, "None", "Adiabatic heating: 0 off; 1 on"),
        iuse_sticky_correction = ParameterInt(
            0, "None", 
            "Sticky temperature correction (reset to top boundary temperature): 0 off; 1 on"
            ),

        # Heat equation - Thermal cond parameters
        itype_conductivity = ParameterInt(
            0, "None", "Integer flag for thermal conductivity model"),
        stype_conductivity = ParameterStr(
            "None", "None", "Thermal conductivity model option name"),

        )
end
