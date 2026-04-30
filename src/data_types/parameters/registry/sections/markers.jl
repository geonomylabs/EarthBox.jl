function get_markers_parameters()::NamedTuple
    return @params (
        # *****************************
        # Markers Parameters
        # *****************************
        
        # Markers - Recycling parameters
        imantle = ParameterInt(
            8, "None", "Material ID of the mantle. This is no longer used and will be removed in the future"),
        # Markers - Distribution parameters (requires constructor args - using defaults),
        iuse_random = ParameterInt(
            0, "None", "Randomize initial marker distribution: 0 off; 1 on"),
        dx_marker = ParameterFloat(
            1000.0, "m", "Initial average marker spacing in horizontal direction in meters"),
        dy_marker = ParameterFloat(
            1000.0, "m", "Initial average marker spacing in vertical direction in meters"),
        nmarkers_cell_x = ParameterFloat(
            3.0, "None", "Number of markers per cell in horizontal direction"),
        nmarkers_cell_y = ParameterFloat(
            3.0, "None", "Number of markers per cell in vertical direction"),
        mxnum = ParameterInt(
            300, "None", "Total number of markers in horizontal direction"),
        mynum = ParameterInt(
            300, "None", "Total number of markers in vertical direction"),
        marknum = ParameterInt(
            90000, "None", "Number of markers"),
        mxstep = ParameterFloat(
            333.33, "m", "Distance between markers in horizontal direction in meters"),
        mystep = ParameterFloat(
            333.33, "m", "Distance between markers in vertical direction in meters"),

        # Markers - Advection parameters
        marker_cell_displ_max = ParameterFloat(
            0.5, "fraction", 
            "Maximum marker displacement as a fraction of cell size in meters"
            ),
        itype_move_markers = ParameterInt(
            4, "None", 
            "Displacement Options: 0 no motion 1 simp. advection; 4 4-th order Runge-Kutta"
            ),
        stype_move_markers = ParameterStr(
            "None", "None", "Displacement options name"),
        iuse_local_adaptive_time_stepping = ParameterInt(
            0, "None",
            "Activate local adaptive time stepping where displacement is calculated at staggered grid nodes. "
            * "The default value is 0 (off). If deactivated then adaptive time stepping using average grid spacing. "
            * "Note: when this flag is 1 the time step can vary significantly between steps, so "
            * "`iuse_fixed_output_counter` is automatically forced to 0 during time-loop initialization "
            * "so that output is triggered by total model time rather than by a fixed step counter."
            ),

        # Markers - Subgrid diffusion parameters
        subgrid_diff_coef_stress = ParameterFloat(
            1.0, "None", "Numerical subgrid stress diffusion coefficient"),
        subgrid_diff_coef_temp = ParameterFloat(
            1.0, "None", "Numerical subgrid temperature diffusion coefficient"),

        )
end
