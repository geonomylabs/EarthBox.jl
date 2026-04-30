function get_carbonate_parameters()::NamedTuple
    return @params (
        iuse_carb = ParameterInt(
            0, "None", "0 = uniform deposition; 1 = deposition on highs"),
        photic_thick_m = ParameterFloat(
            0.0, "meters", "Depth of photic zone"),
        carb_growth_rad = ParameterFloat(
            0.0, "meters", "Carbonate growth radius"),
        carb_base_rate = ParameterFloat(
            0.0, "fraction", "Probability of nucleation"),
        carb_time_myr = ParameterFloat(
            0.0, "Myr", "Carbonate growth delay time relative to model start"),
        carb_jump_time_myr = ParameterFloat(
            0.0, "Myr", "Time of jump in carbonate base rate"),
        carb_base_rate_jump = ParameterFloat(
            0.0, "fraction", "Carbonate base rate jump"),
    )
end
