function get_stokes_continuity_arrays()::NamedTuple
    return @arrays (
        #*************
        # Stokes Continuity Arrays
        #*************

        # Basic Grid Velocity Arrays
        vxb = ArrayData(
            "m/s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : x-component of velocity interpolated to basic grid in m/s.",
        ),
        vyb = ArrayData(
            "m/s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : y-component of velocity interpolated to basic grid in m/s.",
        ),

        # Density Arrays
        rho0 = ArrayData(
            "kg/m^3", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of transport array with density on basic grid in kg/m^3.",
        ),
        rho1 = ArrayData(
            "kg/m^3", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Transport array with density on basic grid in kg/m^3.",
        ),
        rho0_vy = ArrayData(
            "kg/m^3", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Backup of transport array with density on vy grid in kg/m^3.",
        ),
        rho1_vy = ArrayData(
            "kg/m^3", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Transport array with density on vy grid in kg/m^3.",
        ),

        # Residual Stokes Arrays
        resx1 = ArrayData(
            "m/s", ScalarArray2DState, "vx",
            "`(ynum + 1, xnum)` : x-Stokes equation residual on staggered vx grid in m/s.",
        ),
        resnlx = ArrayData(
            "m/s", ScalarArray2DState, "vx",
            "`(ynum + 1, xnum)` : Non-linear x-Stokes equation residual on staggered vx grid in m/s.",
        ),
        resy1 = ArrayData(
            "m/s", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : y-Stokes equation residual on staggered vy grid in m/s.",
        ),
        resnly = ArrayData(
            "m/s", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Non-linear y-Stokes equation residual on staggered vy grid in m/s.",
        ),
        resc1 = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Continuity equation residual on pressure grid in Pa.",
        ),
        resnlc = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Non-linear continuity equation residual on pressure grid in Pa.",
        ),
        resnl = ArrayData(
            "Mixed", SolutionArray1DState, "one-dimensional",
            "Non-linear stokes-continuity system residual.",
        ),

        # Plastic Deformation Arrays
        plastics = ArrayData(
            "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Plastic deformation indicator on basic grid.",
        ),
        plastics0 = ArrayData(
            "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of plastic deformation indicator on basic grid.",
        ),
        plasticn = ArrayData(
            "None", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Plastic deformation indicator on pressure grid.",
        ),
        plasticn0 = ArrayData(
            "None", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Backup of plastic deformation indicator on pressure grid.",
        ),
        plastic_yield = ArrayData(
            "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Plastic yielding marker on basic grid.",
        ),
        yield_error = ArrayData(
            "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Plastic yielding stress error on basic grid.",
        ),
        cohesion_grid = ArrayData(
            "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Cohesion on basic grid interpolated from markers in Pa.",
        ),
        cohesion_grid0 = ArrayData(
            "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of cohesion on basic grid interpolated from markers in Pa.",
        ),
        fric_degrees_grid = ArrayData(
            "Degrees", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Friction angle in degrees on basic grid interpolated from markers in Degrees.",
        ),
        fric_degrees_grid0 = ArrayData(
            "Degrees", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of friction angle in degrees on basic grid in Degrees.",
        ),
        dilatation_grid = ArrayData(
            "Degrees", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Dilatation angle on pressure grid interpolated from markers in Degrees.",
        ),
        dilatation_grid0 = ArrayData(
            "None", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Backup of dilatation angle on pressure grid interpolated from markers.",
        ),
        extractable_meltfrac_grid = ArrayData(
            "Fraction", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Extractable melt fraction on pressure grid in fraction used for melt compaction.",
        ),
        extractable_meltfrac_grid0 = ArrayData(
            "Fraction", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Backup of extractable melt fraction on pressure grid in fraction used for melt compaction.",
        ),

        # RHS Stokes Arrays
        RX1 = ArrayData(
            "m/s or Pa", ScalarArray2DState, "vx",
            "`(ynum + 1, xnum)` : Right-hand side values on vx grid for x-Stokes equation.",
        ),
        RY1 = ArrayData(
            "m/s or Pa", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Right-hand side values on vy grid for y-Stokes equation.",
        ),
        RC1 = ArrayData(
            "", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Right-hand side values on pressure grid for continuity equation.",
        ),
        RHS = ArrayData(
            "", RhsStokesArray1DState, "one-dimensional",
            "Discretized right-hand side array for Stokes-continuity equation",
        ),

        # Shear Modulus Arrays
        mus0 = ArrayData(
            "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of shear modulus for shear stress on basic grid in Pa.",
        ),
        mus1 = ArrayData(
            "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Shear elastic modulus interpolated to basic grid from markers in Pa.",
        ),
        mun0 = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Backup of shear modulus for normal stress on pressure grid in Pa.",
        ),
        mun1 = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Normal elastic modulus interpolated to pressure grid from markers in Pa.",
        ),

        # Pressure Arrays
        pr1_old = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Old pressure solution on pressure grid in Pa.",
        ),
        pr1 = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Pressure solution on pressure grid in Pa.",
        ),

        # Staggered Grid Velocity Arrays
        vy1_old = ArrayData(
            "m/s", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Old y-component of velocity on staggered vy grid in m/s.",
        ),
        vy1 = ArrayData(
            "m/s", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Y-component of velocity on staggered vy grid in m/s.",
        ),
        vx1_old = ArrayData(
            "m/s", ScalarArray2DState, "vx",
            "`(ynum + 1, xnum)` : Old x-component of velocity on staggered vx grid in m/s.",
        ),
        vx1 = ArrayData(
            "m/s", ScalarArray2DState, "vx",
            "`(ynum + 1, xnum)` : X-component of velocity on staggered vx grid in m/s.",
        ),

        # Stokes Solution Arrays
        soluv1_old = ArrayData(
            "m/s", SolutionArray1DState, "one-dimensional",
            "Old Stokes-continuity solution array with velocity and pressure in m/s.",
        ),
        soluv1 = ArrayData(
            "m/s or Pa", SolutionArray1DState, "one-dimensional",
            "Stokes-continuity solution array with velocity and pressure in m/s or Pa.",
        ),

        # Strain Rate and Spin Arrays
        exy = ArrayData(
            "1/s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Shear strain rate in 1/s.",
        ),
        exx = ArrayData(
            "1/s", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Normal strain rate in 1/s.",
        ),
        eii = ArrayData(
            "1/s", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Second invariant of deviatoric strain rate in 1/s.",
        ),
        esp = ArrayData(
            "1/s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Spin in 1/s.",
        ),
        eii_plastic_basic = ArrayData(
            "1/s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Second invariant of plastic strain rate on basic grid in 1/s.",
        ),
        eii_plastic_pressure = ArrayData(
            "1/s", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Second invariant of plastic strain rate on pressure grid in 1/s.",
        ),

        # Velocity Solution Arrays
        vxy_old = ArrayData(
            "m/s", SolutionArray1DState, "one-dimensional",
            "Old solution array with vx and vy solutions in m/s.",
        ),
        vxy = ArrayData(
            "m/s", SolutionArray1DState, "one-dimensional",
            "Solution array with vx and vy solutions in m/s.",
        ),

        # Viscosity Arrays
        etas0 = ArrayData(
            "Pa.s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Old viscoplastic shear viscosity array in Pa.s used for: (1) backup of "
            * "viscoplastic shear viscosity (etas1) interpolated from markers used for "
            * "problematic interpolation and (2) visco-elasto-plastic viscosity coefficients "
            * "used in visco-elastic Stokes equations.",
        ),
        etas1 = ArrayData(
            "Pa.s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Viscoplastic shear viscosity in Pa.s interpolated to basic grid from markers.",
        ),
        etan0 = ArrayData(
            "Pa.s", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Old viscoplastic normal viscosity array in Pa.s used for: (1) "
            * "backup of viscoplastic normal viscosity (etan1) interpolated from markers used "
            * "for problematic interpolation and (2) visco-elasto-plastic viscosity coefficients "
            * "used in visco-elastic Stokes equations.",
        ),
        etan1 = ArrayData(
            "Pa.s", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Viscoplastic normal viscosity in Pa.s interpolated to pressure grid from markers.",
        ),
        eta_flow = ArrayData(
            "Pa.s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Flow viscosity in Pa.s interpolated to basic grid from markers.",
        ),
        eta_flow0 = ArrayData(
            "Pa.s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of flow viscosity on basic grid in Pa.s used for problematic interpolation.",
        ),

        # Stress Change Arrays
        dsxy = ArrayData(
            "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Dual purpose array: First stores post-Stokes viscoelastic stress "
            * "change (sxy2-sxy1) and then updated to final (remaining) shear stress change by "
            * "subtracting relaxed grid nodal-marker stress differences (dsxyn) in Pa.",
        ),
        dsxx = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Dual purpose array: First stores post-Stokes normal "
            * "viscoelastic stress change (sxx2-sxx1) and then updated to final (remaining) "
            * "normal stress change by subtracting relaxed grid nodal-marker stress differences (dsxxn) in Pa.",
        ),
        dsxyn = ArrayData(
            "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Subgrid nodal-marker shear stress change on basic grid interpolated from markers in Pa.",
        ),
        dsxxn = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Subgrid nodal-marker normal stress change on basic grid "
            * "interpolated from markers in Pa.",
        ),

        # Stress Arrays
        sxy0 = ArrayData(
            "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Old deviatoric shear stress array on basic grid in Pa used for: (1) "
            * "storing copy of visco-elasto-plastic shear stress transport array (sxy1) and "
            * "(2) storing product of old visco-elasto-plastic shear stress transport array "
            * "(sxy1) and viscoelastic factor used in RHS terms of visco-elasto-plastic Stokes equations.",
        ),
        sxy1 = ArrayData(
            "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Visco-elasto-plastic deviatoric shear stress in Pa interpolated to basic grid from advected markers.",
        ),
        sxy2 = ArrayData(
            "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Visco-elasto-plastic deviatoric shear stress in Pa on basic grid "
            * "calculated using updated deviatoric strain rates from Stokes velocity solution.",
        ),
        sxx0 = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Old deviatoric normal stress array on basic grid in Pa used for: (1) "
            * "storing copy of visco-elasto-plastic normal stress transport array (sxx1) and "
            * "(2) storing product of old visco-elasto-plastic normal stress transport array "
            * "(sxx1) and viscoelastic factor used in RHS terms of visco-elasto-plastic Stokes equations.",
        ),
        sxx1 = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Visco-elasto-plastic deviatoric normal stress "
            * "interpolated to basic grid from advected markers in Pa.",
        ),
        sxx2 = ArrayData(
            "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Visco-elasto-plastic deviatoric normal stress on "
            * "pressure grid calculated using updated deviatoric strain rates from Stokes velocity solution in Pa.",
        ),
    )
end
