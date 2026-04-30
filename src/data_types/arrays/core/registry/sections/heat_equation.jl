function get_heat_equation_arrays()::NamedTuple
    return @arrays (
        #*************
        # Heat Equation Arrays
        #*************

        # Adiabatic Production Arrays
        ha0 = ArrayData(
            "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of adiabatic heat production term on basic grid.",
        ),
        ha1 = ArrayData(
            "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Adiabatic heat production term (expansivity x temperature) "
            * "interpolated from markers to basic grid.",
        ),

        # Heat Residual Arrays
        rest = ArrayData(
            "", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Residual for heat equation.",
        ),

        # Radiogenic Production Arrays
        hr0 = ArrayData(
            "W/m^3", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of radiogenic heat production on basic grid in W/m^3.",
        ),
        hr1 = ArrayData(
            "W/m^3", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Radiogenic heat production interpolated from markers to "
            * "basic grid in W/m^3.",
        ),

        # RhoCp Arrays
        rhocp0 = ArrayData(
            "kg/m^3 x J/kg/K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of density times heat capacity on basic grid in kg/m^3 x J/kg/K.",
        ),
        rhocp1 = ArrayData(
            "kg/m^3 x J/kg/K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Density times heat capacity interpolated from markers to "
            * "basic grid in kg/m^3 x J/kg/K.",
        ),

        # RHS Heat Arrays
        RT1 = ArrayData(
            "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Right-hand side values at grid node locations.",
        ),
        RHSheat = ArrayData(
            "K", RhsHeatArray1DState, "one-dimensional",
            "`(ynum*xnum)` : Right-hand side array for heat equation.",
        ),

        # Temperature Arrays
        tk0 = ArrayData(
            "K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Old temperature array in Kelvin used for: (1) defining old "
            * "temperature in grids for problematic marker interpolation after being "
            * "copied from last conductive temperature solution (tk2), (2) defining "
            * "old initial temperature when solving transient heat equation after being "
            * "copied from temperature transport array (tk1), and (3) defining old "
            * "temperature from previous iteration of adaptive time stepping heat loop "
            * "after being copied from current conductive temperature solution (tk2).",
        ),
        tk1 = ArrayData(
            "K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Transport array with temperature in Kelvin interpolated from "
            * "advected markers to basic grid used to define initial old temperature "
            * "when solving conductive heat equation.",
        ),
        tk2 = ArrayData(
            "K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : New conductive temperature solution array in Kelvin obtained by "
            * "solving heat conduction equation on basic grid using adaptive time "
            * "stepping heat loop and tk1 as initial condition.",
        ),
        dtk1 = ArrayData(
            "K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Temperature change array in Kelvin used for: (1) defining "
            * "temperature change (tk2 - tk0) on basic grid for iteration of "
            * "adaptive time stepping heat loop, (2) defining temperature change on "
            * "basic grid between new conductive solution and old advective transport "
            * "array temperature (tk2 - tk1). Updated for subgrid diffusion by "
            * "subtracting dtkn.",
        ),
        dtkn = ArrayData(
            "K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Subgrid diffusive temperature change array in Kelvin "
            * "interpolated from markers to basic grid used to update dkt1 for "
            * "subgrid diffusion.",
        ),

        # Thermal Conductivity Arrays
        kt0 = ArrayData(
            "W/m/K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of thermal conductivity on basic grid in W/m/K.",
        ),
        kt1 = ArrayData(
            "W/m/K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Thermal conductivity interpolated from markers to "
            * "basic grid in W/m/K.",
        ),
    )
end
