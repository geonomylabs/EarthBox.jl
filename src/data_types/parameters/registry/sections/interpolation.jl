function get_interpolation_parameters()::NamedTuple
    return @params (
        # *****************************
        # Interpolation Parameters
        # *****************************
        
        # Interpolation parameters
        iuse_initial_temp_interp = ParameterInt(
            0, "None", 
            "Interpolate nodal temperatures back to markers at start: 0 off; 1 on"
            ),
        iuse_harmonic_avg_normal_viscosity = ParameterInt(
            0, "None", 
            "Harmonic averaging of shear viscosity for normal viscosity: 0 off; 1 on"
            ),
        )
end
