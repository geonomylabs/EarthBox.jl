function get_benchmark_parameters()::NamedTuple
    return @params (
        iuse_ramberg_post_processing = ParameterInt(
            0, "None", 
            "Calculate quantities for Rayleigh-Taylor benchmark: 0 = off; 1 = on"),
        iuse_viscous_block_processing = ParameterInt(
            0, "None", 
            "Calculate quantities for viscous block benchmark: 0 = off; 1 = on"),
        iuse_conbox_post_processing = ParameterInt(
            0, "None", 
            "Post-processing for convection in a box: 0 = off; 1 = on"),
    )
end
