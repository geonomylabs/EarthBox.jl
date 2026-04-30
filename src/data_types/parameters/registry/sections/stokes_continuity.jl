function get_stokes_continuity_parameters()::NamedTuple
    return @params (
        # *****************************
        # Stokes-Continuity Parameters
        # *****************************
        
        # Stokes-continuity - Density parameters
        itype_density = ParameterInt(
            0, "None", "Integer ID of density model"),
        stype_density = ParameterStr(
            "Liao14", "None", "String name of density model"),

        # Stokes-continuity - Build parameters (requires constructor args - using defaults),
        ibuild_stokes = ParameterInt(
            1, "None", "0 define full system, 1 only define non-zero elements"),
        hshift_to_vxR = ParameterInt(
            300, "None", 
            "Horizontal shift index for system building equal to (ynum-1)*3"
            ),
        Nstokes = ParameterInt(
            30000, "None", 
            "Number of rows/columns in NxN matrix equal to (xnum-1)*(ynum-1)*3"
            ),
        nonzero_max_stokes = ParameterInt(
            316231, "None", 
            "Maximum number of non-zero values allowed set equal to xnum*ynum*31"
            ),
        pscale = ParameterFloat(
            1.0, "None", "Coefficient for scaling pressure"),
        iuse_interface_stabilization = ParameterInt(
            0, "None", 
            "Flag for interface stabilization: 0 off; 1 on"
            ),

        # Stokes-continuity - Solution norms parameters
        dsoluv1_abs_inf = ParameterFloat(
            1e38, "None", "Infinity norm of solution change"),
        dsoluv1_rel_inf = ParameterFloat(
            1e38, "None", "Relative infinity norm of solution change"),
        dsoluv1_abs_L2 = ParameterFloat(
            1e38, "None", "L2 norm of solution change"),
        dsoluv1_rel_L2 = ParameterFloat(
            1e38, "None", "Relative L2 norm of solution change"),
        dvx1_abs_L2 = ParameterFloat(
            1e38, "None", "L2 norm of vx change"),
        dvx1_rel_L2 = ParameterFloat(
            1e38, "None", "Relative L2 norm of vx change"),
        dvy1_abs_L2 = ParameterFloat(
            1e38, "None", "L2 norm of vy change"),
        dvy1_rel_L2 = ParameterFloat(
            1e38, "None", "Relative L2 norm of vy change"),
        dpr1_abs_L2 = ParameterFloat(
            1e38, "None", "L2 norm of pressure change"),
        dpr1_rel_L2 = ParameterFloat(
            1e38, "None", "Relative L2 norm of pressure change"),
        dvxy_abs_inf = ParameterFloat(
            1e38, "None", "Infinity norm of velocity solution"),
        dvxy_rel_inf = ParameterFloat(
            1e38, "None", "Relative infinity norm of velocity solution"),
        dvxy_abs_L2 = ParameterFloat(
            1e38, "None", "L2 norm of velocity solution"),
        dvxy_rel_L2 = ParameterFloat(
            1e38, "None", "Relative L2 norm of velocity solution"),
        global_yield_error = ParameterFloat(
            1e38, "None", "Global yield stress error"),

        # Stokes-continuity - Residual norms parameters
        resnl_L2_ini = ParameterFloat(
            1e38, "None", "Initial L2 norm of Stokes system"),
        resnl_L2 = ParameterFloat(
            1e38, "None", "Current L2 norm of Stokes system"),
        resnl_rel_L2 = ParameterFloat(
            1e38, "None", "Relative L2 norm of Stokes system"),
        resx_L2 = ParameterFloat(
            1e38, "None", "L2 norm of vx residual"),
        resy_L2 = ParameterFloat(
            1e38, "None", "L2 norm of vy residual"),
        resc_L2 = ParameterFloat(
            1e38, "None", "L2 norm of pressure residual"),
        resx_L2_ini = ParameterFloat(
            1e38, "None", "Initial L2 norm of vx residual"),
        resy_L2_ini = ParameterFloat(
            1e38, "None", "Initial L2 norm of vy residual"),
        resc_L2_ini = ParameterFloat(
            1e38, "None", "Initial L2 norm of pressure residual"),
        resnlx_L2 = ParameterFloat(
            1e38, "None", "Non-linear L2 norm of vx residual"),
        resnly_L2 = ParameterFloat(
            1e38, "None", "Non-linear L2 norm of vy residual"),
        resnlc_L2 = ParameterFloat(
            1e38, "None", "Non-linear L2 norm of pressure residual"),
        resnlx_rel_L2 = ParameterFloat(
            1e38, "None", "Relative non-linear L2 norm of vx residual"),
        resnly_rel_L2 = ParameterFloat(
            1e38, "None", "Relative non-linear L2 norm of vy residual"),
        resnlc_rel_L2 = ParameterFloat(
            1e38, "None", "Relative non-linear L2 norm of pressure residual"),
        resnlx_L2_ini = ParameterFloat(
            1e38, "None", "Initial non-linear L2 norm of vx residual"),
        resnly_L2_ini = ParameterFloat(
            1e38, "None", "Initial non-linear L2 norm of vy residual"),
        resnlc_L2_ini = ParameterFloat(
            1e38, "None", "Initial non-linear L2 norm of pressure residual"),

        # Stokes-continuity - Picard parameters
        nglobal = ParameterInt(
            1, "None", 
            "Maximum number of global plasticity iterations (a.k.a. Picard iterations)"
            ),
        tolerance_picard = ParameterFloat(
            1e-2, "None", 
            "Convergence criterion for the global plasticity iterations (a.k.a. Picard iterations)"
            ),
        iconverge = ParameterInt(
            0, "None", "Convergence flag"),
        iglobal = ParameterInt(
            0, "None", "Iteration counter for global stokes loop"),
        itype_global = ParameterInt(
            1, "None", "Type of global loop: 0 = marker-base; 1 = node-based"),
        stype_global = ParameterStr(
            "None", "None", "Global loop option name"),

        # Stokes-continuity - Velocity calc options parameters
        itype_velocity = ParameterInt(0, "None", 
            "Velocity calculation options: -1 all zeros; 0 solve equations; 1 solid body rotation"),
        stype_velocity = ParameterStr(
            "None", "None", "Velocity calculation option name"),
        # Stokes-continuity - Output parameters
        outtest = ParameterInt(
            0, "None", "0 off; 1 output grid, vis., RHS, and bc as text files"),
        outtest2 = ParameterInt(
            0, "None", "1 output sparse matrix, rhs and solution vector as text files"),
        outtest3 = ParameterInt(
            0, "None", "1 output sparse matrix, rhs and solution vector as npy files"),

        )
end
