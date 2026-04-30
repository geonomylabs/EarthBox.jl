function get_grid_parameters()::NamedTuple
    return @params (
        # *****************************
        # Grid Parameters
        # *****************************

        # Grid - Geometry parameters (requires constructor args - using defaults),
        ynum = ParameterInt(
            101, "None", "Number of basic grid points in y-direction"),
        xnum = ParameterInt(
            101, "None", "Number of basic grid points in x-direction"),
        ysize = ParameterFloat(
            100000.0, "m", "Height of model in meters"),
        xsize = ParameterFloat(
            100000.0, "m", "Width of model in meters"),
        ystpavg = ParameterFloat(
            0.0, "m", "Average spacing of basic grid in y-direction in meters"),
        xstpavg = ParameterFloat(
            0.0, "m", "Average spacing of basic grid in x-direction in meters"),
        ymin = ParameterFloat(
            0.0, "m", "Minimum y-coordinate of model domain in meters"),
        ymax = ParameterFloat(
            100000.0, "m", "Maximum y-coordinate of model domain in meters"),
        xmin = ParameterFloat(
            0.0, "m", "Minimum x-coordinate of model domain in meters"),
        xmax = ParameterFloat(
            100000.0, "m", "Maximum x-coordinate of model domain in meters"),
        xsize_start = ParameterFloat(
            100000.0, "m", "Initial width of model in meters"),

        # Grid - Options parameters
        itype_grid = ParameterInt(
            0, "None", "Grid option ID"),
        stype_grid = ParameterStr(
            "None", "None", "Grid option name"),

        # Grid - Refinement parameters
        dx_highres = ParameterFloat(
            NaN, "m", 
            "Constant horizontal grid resolution in high-resolution area in meters"
            ),
        dx_lowres = ParameterFloat(
            NaN, "m", 
            "Average horizontal grid resolution in low-resolution area in meters"
            ),
        xo_highres = ParameterFloat(
            NaN, "m", "x-location of first node of high-resolution area in meters"),
        ixo_highres = ParameterInt(
            -9999, "None", "x-index of first node of high-resolution area"),
        xf_highres = ParameterFloat(
            NaN, "m", "x-location of last node of high-resolution area in meters"),
        dy_highres = ParameterFloat(
            NaN, "m", "Constant vertical grid resolution in high-resolution area in meters"),
        dy_lowres = ParameterFloat(
            NaN, "m", "Average vertical grid resolution in low-resolution area in meters"),
        yf_highres = ParameterFloat(
            NaN, "m", "y-location of last node of high-resolution area in meters"),
        iuse_trench = ParameterInt(
            0, "None", "Flag to use trench location for refinement"),
        trench_location = ParameterFloat(
            NaN, "m", "x-location of trench in meters"),
        iuse_refinement_delay = ParameterInt(
            0, "None", "Flag to delay grid refinement"),
        refinement_time = ParameterFloat(
            NaN, "Myr", "Time to start grid refinement in Myr"),
        refinement_flag = ParameterInt(
            0, "None", "Flag indicating if grid is currently refined"),
        iuse_refinement_gap = ParameterInt(
            0, "None", "Flag to temporarily disable refinement"),
        refinement_gap_start_time = ParameterFloat(
            NaN, "Myr", "Time to start refinement gap in Myr"),
        refinement_gap_end_time = ParameterFloat(
            NaN, "Myr", "Time to end refinement gap in Myr"),

        )
end
