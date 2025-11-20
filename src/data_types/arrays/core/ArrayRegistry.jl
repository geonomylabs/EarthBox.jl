module ArrayRegistry

import EarthBox.Arrays.ArrayTypes.GridArray1D: GridArray1DState
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.Arrays.ArrayTypes.RhsHeatArray1D: RhsHeatArray1DState
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.Arrays.ArrayTypes.MarkerArrayInt1D: MarkerArrayInt1DState
import EarthBox.Arrays.ArrayTypes.SolutionArray1D: SolutionArray1DState
import EarthBox.Arrays.ArrayTypes.RhsStokesArray1D: RhsStokesArray1DState
import EarthBox.Arrays.ArrayTypes.BcArrayFloat: BcArrayFloatState
import EarthBox.Arrays.ArrayTypes.InternalBcArrayInt: InternalBcArrayIntState
import EarthBox.Arrays.ArrayTypes.InternalBcArrayFloat: InternalBcArrayFloatState
import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray1D
import EarthBox.EarthBoxDtypes: AbstractEarthBoxArray2D

struct ArrayData
    name::String
    units::String
    type::Union{Type{<:AbstractEarthBoxArray1D}, Type{<:AbstractEarthBoxArray2D}}
    grid_type::String
    description::String
end

function get_eb_arrays()::NamedTuple
    return merge_named_tuples(
        get_grid_arrays(),
        get_heat_equation_arrays(),
        get_interpolation_arrays(),
        get_markers_arrays(),
        get_melting_arrays(),
        get_stokes_continuity_arrays(),
        get_bcs_arrays()
    )
end

function merge_named_tuples(tuples...)
    result = tuples[1]
    ntuples = length(tuples)
    for i in 2:ntuples
        result = merge(result, tuples[i])
    end
    return result
end

function get_grid_arrays()::NamedTuple
    return (
        #*************
        # Grids Arrays
        #*************

        # Basic Grid Arrays
        gridx_b = ArrayData(
            "gridx_b", "m", GridArray1DState, "one-dimensional",
            "`(xnum)` : x-coordinates of basic grid nodes in meters.",
        ),
        xstp_b = ArrayData(
            "xstp_b", "m", GridArray1DState, "one-dimensional",
            "`(xnum - 1)` : Width of cells in x-direction for basic grid in meters.",
        ),
        gridy_b = ArrayData(
            "gridy_b", "m", GridArray1DState, "one-dimensional",
            "`(ynum)` : y-coordinates of basic grid nodes in meters.",
        ),
        ystp_b = ArrayData(
            "ystp_b", "m", GridArray1DState, "one-dimensional",
            "`(ynum - 1)` : Width of cells in y-direction for basic grid in meters.",
        ),
        # Pressure Grid Arrays
        gridy_pr = ArrayData(
            "gridy_pr", "m", GridArray1DState, "one-dimensional",
            "`(ynum - 1)` : y-locations of pressure grid nodes in meters.",
        ),
        ystp_pr = ArrayData(
            "ystp_pr", "m", GridArray1DState, "one-dimensional",
            "`(ynum - 2)` : Width of pressure cells in y-direction in meters.",
        ),
        gridx_pr = ArrayData(
            "gridx_pr", "m", GridArray1DState, "one-dimensional",
            "`(xnum - 1)` : x-locations of pressure grid nodes in meters.",
        ),
        xstp_pr = ArrayData(
            "xstp_pr", "m", GridArray1DState, "one-dimensional",
            "`(xnum - 2)` : Width of pressure cells in x-direction in meters.",
        ),
        # Staggered Vx Grid Arrays
        gridy_vx = ArrayData(
            "gridy_vx", "m", GridArray1DState, "one-dimensional",
            "`(ynum + 1)` : y-locations of staggered Vx grid nodes including ghost nodes in meters.",
        ),
        ystp_vx = ArrayData(
            "ystp_vx", "m", GridArray1DState, "one-dimensional",
            "`(ynum)` : Width of staggered vx grid cells in y-direction in meters.",
        ),
        # Staggered Vy Grid Arrays
        gridx_vy = ArrayData(
            "gridx_vy", "m", GridArray1DState, "one-dimensional",
            "`(xnum + 1)` : x-locations of staggered Vy grid nodes including ghost nodes in meters.",
        ),
        xstp_vy = ArrayData(
            "xstp_vy", "m", GridArray1DState, "one-dimensional",
            "`(xnum)` : Width of staggered vy grid cells in x-direction in meters.",
        ),
    )

end

function get_heat_equation_arrays()::NamedTuple
    return (
        #*************
        # Heat Equation Arrays
        #*************

        # Adiabatic Production Arrays
        ha0 = ArrayData(
            "ha0", "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of adiabatic heat production term on basic grid.",
        ),
        ha1 = ArrayData(
            "ha1", "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Adiabatic heat production term (expansivity x temperature) "
            * "interpolated from markers to basic grid.",
        ),

        # Heat Residual Arrays
        rest = ArrayData(
            "rest", "", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Residual for heat equation.",
        ),

        # Radiogenic Production Arrays
        hr0 = ArrayData(
            "hr0", "W/m^3", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of radiogenic heat production on basic grid in W/m^3.",
        ),
        hr1 = ArrayData(
            "hr1", "W/m^3", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Radiogenic heat production interpolated from markers to "
            * "basic grid in W/m^3.",
        ),

        # RhoCp Arrays
        rhocp0 = ArrayData(
            "rhocp0", "kg/m^3 x J/kg/K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of density times heat capacity on basic grid in kg/m^3 x J/kg/K.",
        ),
        rhocp1 = ArrayData(
            "rhocp1", "kg/m^3 x J/kg/K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Density times heat capacity interpolated from markers to "
            * "basic grid in kg/m^3 x J/kg/K.",
        ),

        # RHS Heat Arrays
        RT1 = ArrayData(
            "RT1", "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Right-hand side values at grid node locations.",
        ),
        RHSheat = ArrayData(
            "RHSheat", "K", RhsHeatArray1DState, "one-dimensional",
            "`(ynum*xnum)` : Right-hand side array for heat equation.",
        ),

        # Temperature Arrays
        tk0 = ArrayData(
            "tk0", "K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Old temperature array in Kelvin used for: (1) defining old "
            * "temperature in grids for problematic marker interpolation after being "
            * "copied from last conductive temperature solution (tk2), (2) defining "
            * "old initial temperature when solving transient heat equation after being "
            * "copied from temperature transport array (tk1), and (3) defining old "
            * "temperature from previous iteration of adaptive time stepping heat loop "
            * "after being copied from current conductive temperature solution (tk2).",
        ),
        tk1 = ArrayData(
            "tk1", "K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Transport array with temperature in Kelvin interpolated from "
            * "advected markers to basic grid used to define initial old temperature "
            * "when solving conductive heat equation.",
        ),
        tk2 = ArrayData(
            "tk2", "K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : New conductive temperature solution array in Kelvin obtained by "
            * "solving heat conduction equation on basic grid using adaptive time "
            * "stepping heat loop and tk1 as initial condition.",
        ),
        dtk1 = ArrayData(
            "dtk1", "K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Temperature change array in Kelvin used for: (1) defining "
            * "temperature change (tk2 - tk0) on basic grid for iteration of "
            * "adaptive time stepping heat loop, (2) defining temperature change on "
            * "basic grid between new conductive solution and old advective transport "
            * "array temperature (tk2 - tk1). Updated for subgrid diffusion by "
            * "subtracting dtkn.",
        ),
        dtkn = ArrayData(
            "dtkn", "K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Subgrid diffusive temperature change array in Kelvin "
            * "interpolated from markers to basic grid used to update dkt1 for "
            * "subgrid diffusion.",
        ),

        # Thermal Conductivity Arrays
        kt0 = ArrayData(
            "kt0", "W/m/K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of thermal conductivity on basic grid in W/m/K.",
        ),
        kt1 = ArrayData(
            "kt1", "W/m/K", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Thermal conductivity interpolated from markers to "
            * "basic grid in W/m/K.",
        ),
    )
end


function get_interpolation_arrays()::NamedTuple
    return (
        #*************
        # Interpolation Arrays
        #*************

        # Grid Weights Arrays
        wtnodes = ArrayData(
            "wtnodes", "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Summed marker weight factors on basic grid nodes "
            * "using inclusive method with maximum search radius for markers.",
        ),
        wtetas = ArrayData(
            "wtetas", "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Summed marker weight factors on basic grid nodes "
            * "using exclusive method with shorter search radius for markers.",
        ),
        wtetan = ArrayData(
            "wtetan", "None", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Summed marker weight factors for normal "
            * "viscosity on pressure nodes.",
        ),
        wtnodes_vy = ArrayData(
            "wtnodes_vy", "None", ScalarArray2DState, "vy",
            "`(ynum + 1, xnum)` : Summed marker weight factors on Vy staggered "
            * "grid nodes using inclusive method with maximum search radius.",
        ),

        # Marker Weights Arrays
        marker_wtforULnode = ArrayData(
            "marker_wtforULnode", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for upper-left basic grid node.",
        ),
        marker_wtforLLnode = ArrayData(
            "marker_wtforLLnode", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for lower-left basic grid node.",
        ),
        marker_wtforURnode = ArrayData(
            "marker_wtforURnode", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for upper-right basic grid node.",
        ),
        marker_wtforLRnode = ArrayData(
            "marker_wtforLRnode", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for lower-right basic grid node.",
        ),
        marker_wtforCnode = ArrayData(
            "marker_wtforCnode", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for central pressure grid node.",
        ),
        marker_wtforULnodeVy = ArrayData(
            "marker_wtforULnodeVy", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for upper-left Vy staggered grid node.",
        ),
        marker_wtforLLnodeVy = ArrayData(
            "marker_wtforLLnodeVy", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for lower-left Vy staggered grid node.",
        ),
        marker_wtforURnodeVy = ArrayData(
            "marker_wtforURnodeVy", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for upper-right Vy staggered grid node.",
        ),
        marker_wtforLRnodeVy = ArrayData(
            "marker_wtforLRnodeVy", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for lower-right Vy staggered grid node.",
        ),
    )
end

function get_markers_arrays()::NamedTuple
    return (
        #*************
        # Markers Arrays
        #*************

        # Grid Marker Relationship Arrays
        marker_xn = ArrayData(
            "marker_xn", "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Horizontal index of upper left node of basic grid cell containing marker.",
        ),
        marker_yn = ArrayData(
            "marker_yn", "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Vertical index of upper left node of basic grid cell containing marker.",
        ),
        marker_dx = ArrayData(
            "marker_dx", "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Normalized x-distance between marker and upper-left basic grid node in meters.",
        ),
        marker_dy = ArrayData(
            "marker_dy", "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Normalized y-distance between marker and upper-left basic grid node in meters.",
        ),
        marker_xn_vy = ArrayData(
            "marker_xn_vy", "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Horizontal index of upper left node of staggered vy grid cell containing marker.",
        ),
        marker_yn_vy = ArrayData(
            "marker_yn_vy", "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Vertical index of upper left node of staggered vy grid cell containing marker.",
        ),
        marker_dx_vy = ArrayData(
            "marker_dx_vy", "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Normalized x-distance between marker and upper-left grid node of staggered vy grid in meters.",
        ),
        marker_dy_vy = ArrayData(
            "marker_dy_vy", "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Normalized y-distance between marker and upper-left grid node of staggered vy grid in meters.",
        ),

        # Location Arrays
        marker_x = ArrayData(
            "marker_x", "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : x-coordinate of marker in meters.",
        ),
        marker_y = ArrayData(
            "marker_y", "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : y-coordinate of marker in meters.",
        ),

        # Material Arrays
        marker_matid = ArrayData(
            "marker_matid", "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Material ID of marker.",
        ),
        marker_rho = ArrayData(
            "marker_rho", "kg/m^3", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker density in kg/m^3.",
        ),
        marker_porosity_initial = ArrayData(
            "marker_porosity_initial", "fraction", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Initial porosity at sediment-water interface used in Athy's law in fraction.",
        ),
        marker_decay_depth = ArrayData(
            "marker_decay_depth", "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Porosity decay depth used in Athy's law in meters.",
        ),
        marker_max_burial_depth = ArrayData(
            "marker_max_burial_depth", "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Maximum burial depth of marker in meters.",
        ),
        marker_serpentinization = ArrayData(
            "marker_serpentinization", "fraction", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Serpentinization fraction of marker in fraction.",
        ),
        marker_serpentinization_heat_production = ArrayData(
            "marker_serpentinization_heat_production", "W/m^3", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Serpentinization heat production of marker in W/m^3.",
        ),

        # Melt Arrays
        marker_meltfrac = ArrayData(
            "marker_meltfrac", "fraction", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Melt fraction of marker in fraction.",
        ),
        marker_extracted_meltfrac = ArrayData(
            "marker_extracted_meltfrac", "fraction", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Extracted melt fraction of marker in fraction.",
        ),
        marker_extractable_meltfrac = ArrayData(
            "marker_extractable_meltfrac", "fraction", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Extractable melt fraction of marker in fraction.",
        ),

        # Pressure Arrays
        marker_pr = ArrayData(
            "marker_pr", "Pa", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pressure of marker in Pa.",
        ),

        # Rheology Arrays
        marker_eta = ArrayData(
            "marker_eta", "Pa.s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker viscoplastic viscosity in Pa.s.",
        ),
        marker_fric_ini = ArrayData(
            "marker_fric_ini", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Initial friction coefficient (sine of friction angle) of marker.",
        ),
        marker_fric = ArrayData(
            "marker_fric", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Friction coefficient (sine of friction angle) of marker.",
        ),
        marker_cohesion = ArrayData(
            "marker_cohesion", "Pa", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Cohesion of marker in Pa.",
        ),
        marker_preexp = ArrayData(
            "marker_preexp", "1/s/MPa^n", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-exponential factor for dislocation creep of marker in 1/s/MPa^n.",
        ),
        marker_pfailure = ArrayData(
            "marker_pfailure", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Flag indicating plastic failure of marker.",
        ),
        marker_eta_flow = ArrayData(
            "marker_eta_flow", "Pa.s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker viscosity from flow law in Pa.s.",
        ),
        marker_dilatation_angle = ArrayData(
            "marker_dilatation_angle", "Degrees", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker dilatation angle in Degrees.",
        ),

        # Strain Arrays
        marker_GII = ArrayData(
            "marker_GII", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Accumulated strain of marker.",
        ),
        marker_strain_plastic = ArrayData(
            "marker_strain_plastic", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Accumulated plastic strain of marker.",
        ),
        marker_exx = ArrayData(
            "marker_exx", "1/s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Normal strain rate of marker in 1/s.",
        ),
        marker_exy = ArrayData(
            "marker_exy", "1/s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Shear strain rate of marker in 1/s.",
        ),
        marker_sr_ratio = ArrayData(
            "marker_sr_ratio", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Ratio of nodal-stress-based to nodal-velocity-based strain rate invariant for each marker.",
        ),
        marker_strain_rate_plastic = ArrayData(
            "marker_strain_rate_plastic", "1/s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Plastic strain rate of marker in 1/s.",
        ),
        marker_melt_damage = ArrayData(
            "marker_melt_damage", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Damage factor of marker associated with melt intrusion.",
        ),

        # Stratigraphy Arrays
        marker_age = ArrayData(
            "marker_age", "Myr", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Age of marker in million years.",
        ),

        # Stress Arrays
        marker_sxx = ArrayData(
            "marker_sxx", "Pa", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Deviatoric normal stress of marker in Pa.",
        ),
        marker_sxy = ArrayData(
            "marker_sxy", "Pa", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Deviatoric shear stress of marker in Pa.",
        ),

        # Thermal Arrays
        marker_TK = ArrayData(
            "marker_TK", "K", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Temperature of marker in Kelvin.",
        ),
        marker_rhocp = ArrayData(
            "marker_rhocp", "J/K/m^3", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker density multiplied by heat capacity in J/K/m^3.",
        ),
        marker_kt = ArrayData(
            "marker_kt", "W/m/K", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker thermal conductivity in W/m/K.",
        ),
        marker_ha = ArrayData(
            "marker_ha", "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker adiabatic heating term (expansivity x temperature).",
        ),
    )
end

function get_melting_arrays()::NamedTuple
    return (
        #*************
        # Melting Arrays
        #*************

        # Extraction Arrays
        xstart_drainage = ArrayData(
            "xstart_drainage", "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Current x-locations in meters of the start of the drainage area. "
            * "The drainage area is the region where the melt migrates to the "
            * "shallowest position.",
        ),
        xend_drainage = ArrayData(
            "xend_drainage", "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Current x-locations in meters of the end of the drainage area. "
            * "The drainage area is the region where the melt migrates to the "
            * "shallowest position.",
        ),
        melt_residuals = ArrayData(
            "melt_residuals", "m^3", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Residual melt volume in m^3 remaining in the extrusion volume "
            * "after the melt has been extracted in drainage basin.",
        ),
        extrusion_volumes = ArrayData(
            "extrusion_volumes", "m^3", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Volume of extracted melt in m^3 that will be extruded to surface.",
        ),
        xmid_molten_zones = ArrayData(
            "xmid_molten_zones", "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : X-locations in meters of the midpoints of the molten region.",
        ),
        ytop_molten_zones = ArrayData(
            "ytop_molten_zones", "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Y-locations in meters of the top of the molten region.",
        ),
        width_molten_zones = ArrayData(
            "width_molten_zones", "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Width of the molten region in drainage basin in meters.",
        ),
        avg_shallow_partial_melt_xcoors = ArrayData(
            "avg_shallow_partial_melt_xcoors", "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Average x-coordinates of the shallowest partial melt "
            * "marker in drainage basins in meters.",
        ),
        avg_shallow_partial_melt_ycoors = ArrayData(
            "avg_shallow_partial_melt_ycoors", "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Average y-coordinates of the shallowest partial melt "
            * "marker in drainage basins in meters.",
        ),
        xstart_drainage_o = ArrayData(
            "xstart_drainage_o", "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Old x-locations of the start of the drainage area. "
            * "The drainage area is the region where the melt migrates to the "
            * "shallowest position in meters.",
        ),
        xend_drainage_o = ArrayData(
            "xend_drainage_o", "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Old x-locations of the end of the drainage area. "
            * "The drainage area is the region where the melt migrates to the "
            * "shallowest position in meters.",
        ),
    )
end

function get_stokes_continuity_arrays()::NamedTuple
    return (
        #*************
        # Stokes Continuity Arrays
        #*************

        # Basic Grid Velocity Arrays
        vxb = ArrayData(
            "vxb", "m/s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : x-component of velocity interpolated to basic grid in m/s.",
        ),
        vyb = ArrayData(
            "vyb", "m/s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : y-component of velocity interpolated to basic grid in m/s.",
        ),

        # Density Arrays
        rho0 = ArrayData(
            "rho0", "kg/m^3", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of transport array with density on basic grid in kg/m^3.",
        ),
        rho1 = ArrayData(
            "rho1", "kg/m^3", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Transport array with density on basic grid in kg/m^3.",
        ),
        rho0_vy = ArrayData(
            "rho0_vy", "kg/m^3", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Backup of transport array with density on vy grid in kg/m^3.",
        ),
        rho1_vy = ArrayData(
            "rho1_vy", "kg/m^3", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Transport array with density on vy grid in kg/m^3.",
        ),

        # Residual Stokes Arrays
        resx1 = ArrayData(
            "resx1", "m/s", ScalarArray2DState, "vx",
            "`(ynum + 1, xnum)` : x-Stokes equation residual on staggered vx grid in m/s.",
        ),
        resnlx = ArrayData(
            "resnlx", "m/s", ScalarArray2DState, "vx",
            "`(ynum + 1, xnum)` : Non-linear x-Stokes equation residual on staggered vx grid in m/s.",
        ),
        resy1 = ArrayData(
            "resy1", "m/s", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : y-Stokes equation residual on staggered vy grid in m/s.",
        ),
        resnly = ArrayData(
            "resnly", "m/s", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Non-linear y-Stokes equation residual on staggered vy grid in m/s.",
        ),
        resc1 = ArrayData(
            "resc1", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Continuity equation residual on pressure grid in Pa.",
        ),
        resnlc = ArrayData(
            "resnlc", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Non-linear continuity equation residual on pressure grid in Pa.",
        ),
        resnl = ArrayData(
            "resnl", "Mixed", SolutionArray1DState, "one-dimensional",
            "Non-linear stokes-continuity system residual.",
        ),

        # Plastic Deformation Arrays
        plastics = ArrayData(
            "plastics", "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Plastic deformation indicator on basic grid.",
        ),
        plastics0 = ArrayData(
            "plastics0", "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of plastic deformation indicator on basic grid.",
        ),
        plasticn = ArrayData(
            "plasticn", "None", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Plastic deformation indicator on pressure grid.",
        ),
        plasticn0 = ArrayData(
            "plasticn0", "None", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Backup of plastic deformation indicator on pressure grid.",
        ),
        plastic_yield = ArrayData(
            "plastic_yield", "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Plastic yielding marker on basic grid.",
        ),
        yield_error = ArrayData(
            "yield_error", "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Plastic yielding stress error on basic grid.",
        ),
        cohesion_grid = ArrayData(
            "cohesion_grid", "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Cohesion on basic grid interpolated from markers in Pa.",
        ),
        cohesion_grid0 = ArrayData(
            "cohesion_grid0", "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of cohesion on basic grid interpolated from markers in Pa.",
        ),
        fric_degrees_grid = ArrayData(
            "fric_degrees_grid", "Degrees", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Friction angle in degrees on basic grid interpolated from markers in Degrees.",
        ),
        fric_degrees_grid0 = ArrayData(
            "fric_degrees_grid0", "Degrees", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of friction angle in degrees on basic grid in Degrees.",
        ),
        dilatation_grid = ArrayData(
            "dilatation_grid", "Degrees", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Dilatation angle on pressure grid interpolated from markers in Degrees.",
        ),
        dilatation_grid0 = ArrayData(
            "dilatation_grid0", "None", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Backup of dilatation angle on pressure grid interpolated from markers.",
        ),

        # RHS Stokes Arrays
        RX1 = ArrayData(
            "RX1", "m/s or Pa", ScalarArray2DState, "vx",
            "`(ynum + 1, xnum)` : Right-hand side values on vx grid for x-Stokes equation.",
        ),
        RY1 = ArrayData(
            "RY1", "m/s or Pa", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Right-hand side values on vy grid for y-Stokes equation.",
        ),
        RC1 = ArrayData(
            "RC1", "", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Right-hand side values on pressure grid for continuity equation.",
        ),
        RHS = ArrayData(
            "RHS", "", RhsStokesArray1DState, "one-dimensional",
            "Discretized right-hand side array for Stokes-continuity equation",
        ),

        # Shear Modulus Arrays
        mus0 = ArrayData(
            "mus0", "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of shear modulus for shear stress on basic grid in Pa.",
        ),
        mus1 = ArrayData(
            "mus1", "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Shear elastic modulus interpolated to basic grid from markers in Pa.",
        ),
        mun0 = ArrayData(
            "mun0", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Backup of shear modulus for normal stress on pressure grid in Pa.",
        ),
        mun1 = ArrayData(
            "mun1", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Normal elastic modulus interpolated to pressure grid from markers in Pa.",
        ),

        # Pressure Arrays
        pr1_old = ArrayData(
            "pr1_old", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Old pressure solution on pressure grid in Pa.",
        ),
        pr1 = ArrayData(
            "pr1", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Pressure solution on pressure grid in Pa.",
        ),

        # Staggered Grid Velocity Arrays
        vy1_old = ArrayData(
            "vy1_old", "m/s", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Old y-component of velocity on staggered vy grid in m/s.",
        ),
        vy1 = ArrayData(
            "vy1", "m/s", ScalarArray2DState, "vy",
            "`(ynum, xnum + 1)` : Y-component of velocity on staggered vy grid in m/s.",
        ),
        vx1_old = ArrayData(
            "vx1_old", "m/s", ScalarArray2DState, "vx",
            "`(ynum + 1, xnum)` : Old x-component of velocity on staggered vx grid in m/s.",
        ),
        vx1 = ArrayData(
            "vx1", "m/s", ScalarArray2DState, "vx",
            "`(ynum + 1, xnum)` : X-component of velocity on staggered vx grid in m/s.",
        ),

        # Stokes Solution Arrays
        soluv1_old = ArrayData(
            "soluv1_old", "m/s", SolutionArray1DState, "one-dimensional",
            "Old Stokes-continuity solution array with velocity and pressure in m/s.",
        ),
        soluv1 = ArrayData(
            "soluv1", "m/s or Pa", SolutionArray1DState, "one-dimensional",
            "Stokes-continuity solution array with velocity and pressure in m/s or Pa.",
        ),

        # Strain Rate and Spin Arrays
        exy = ArrayData(
            "exy", "1/s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Shear strain rate in 1/s.",
        ),
        exx = ArrayData(
            "exx", "1/s", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Normal strain rate in 1/s.",
        ),
        eii = ArrayData(
            "eii", "1/s", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Second invariant of deviatoric strain rate in 1/s.",
        ),
        esp = ArrayData(
            "esp", "1/s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Spin in 1/s.",
        ),
        eii_plastic_basic = ArrayData(
            "eii_plastic_basic", "1/s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Second invariant of plastic strain rate on basic grid in 1/s.",
        ),
        eii_plastic_pressure = ArrayData(
            "eii_plastic_pressure", "1/s", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Second invariant of plastic strain rate on pressure grid in 1/s.",
        ),

        # Velocity Solution Arrays
        vxy_old = ArrayData(
            "vxy_old", "m/s", SolutionArray1DState, "one-dimensional",
            "Old solution array with vx and vy solutions in m/s.",
        ),
        vxy = ArrayData(
            "vxy", "m/s", SolutionArray1DState, "one-dimensional",
            "Solution array with vx and vy solutions in m/s.",
        ),

        # Viscosity Arrays
        etas0 = ArrayData(
            "etas0", "Pa.s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Old viscoplastic shear viscosity array in Pa.s used for: (1) backup of "
            * "viscoplastic shear viscosity (etas1) interpolated from markers used for "
            * "problematic interpolation and (2) visco-elasto-plastic viscosity coefficients "
            * "used in visco-elastic Stokes equations.",
        ),
        etas1 = ArrayData(
            "etas1", "Pa.s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Viscoplastic shear viscosity in Pa.s interpolated to basic grid from markers.",
        ),
        etan0 = ArrayData(
            "etan0", "Pa.s", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Old viscoplastic normal viscosity array in Pa.s used for: (1) "
            * "backup of viscoplastic normal viscosity (etan1) interpolated from markers used "
            * "for problematic interpolation and (2) visco-elasto-plastic viscosity coefficients "
            * "used in visco-elastic Stokes equations.",
        ),
        etan1 = ArrayData(
            "etan1", "Pa.s", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Viscoplastic normal viscosity in Pa.s interpolated to pressure grid from markers.",
        ),
        eta_flow = ArrayData(
            "eta_flow", "Pa.s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Flow viscosity in Pa.s interpolated to basic grid from markers.",
        ),
        eta_flow0 = ArrayData(
            "eta_flow0", "Pa.s", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Backup of flow viscosity on basic grid in Pa.s used for problematic interpolation.",
        ),

        # Stress Change Arrays
        dsxy = ArrayData(
            "dsxy", "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Dual purpose array: First stores post-Stokes viscoelastic stress "
            * "change (sxy2-sxy1) and then updated to final (remaining) shear stress change by "
            * "subtracting relaxed grid nodal-marker stress differences (dsxyn) in Pa.",
        ),
        dsxx = ArrayData(
            "dsxx", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Dual purpose array: First stores post-Stokes normal "
            * "viscoelastic stress change (sxx2-sxx1) and then updated to final (remaining) "
            * "normal stress change by subtracting relaxed grid nodal-marker stress differences (dsxxn) in Pa.",
        ),
        dsxyn = ArrayData(
            "dsxyn", "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Subgrid nodal-marker shear stress change on basic grid interpolated from markers in Pa.",
        ),
        dsxxn = ArrayData(
            "dsxxn", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Subgrid nodal-marker normal stress change on basic grid "
            * "interpolated from markers in Pa.",
        ),

        # Stress Arrays
        sxy0 = ArrayData(
            "sxy0", "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Old deviatoric shear stress array on basic grid in Pa used for: (1) "
            * "storing copy of visco-elasto-plastic shear stress transport array (sxy1) and "
            * "(2) storing product of old visco-elasto-plastic shear stress transport array "
            * "(sxy1) and viscoelastic factor used in RHS terms of visco-elasto-plastic Stokes equations.",
        ),
        sxy1 = ArrayData(
            "sxy1", "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Visco-elasto-plastic deviatoric shear stress in Pa interpolated to basic grid from advected markers.",
        ),
        sxy2 = ArrayData(
            "sxy2", "Pa", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Visco-elasto-plastic deviatoric shear stress in Pa on basic grid "
            * "calculated using updated deviatoric strain rates from Stokes velocity solution.",
        ),
        sxx0 = ArrayData(
            "sxx0", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Old deviatoric normal stress array on basic grid in Pa used for: (1) "
            * "storing copy of visco-elasto-plastic normal stress transport array (sxx1) and "
            * "(2) storing product of old visco-elasto-plastic normal stress transport array "
            * "(sxx1) and viscoelastic factor used in RHS terms of visco-elasto-plastic Stokes equations.",
        ),
        sxx1 = ArrayData(
            "sxx1", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Visco-elasto-plastic deviatoric normal stress "
            * "interpolated to basic grid from advected markers in Pa.",
        ),
        sxx2 = ArrayData(
            "sxx2", "Pa", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Visco-elasto-plastic deviatoric normal stress on "
            * "pressure grid calculated using updated deviatoric strain rates from Stokes velocity solution in Pa.",
        ),
    )
end

function get_bcs_arrays()::NamedTuple
    return (
        #*************
        # Boundary Conditions Arrays
        #*************

        # Temperature Boundary Conditions Arrays
        btopt = ArrayData(
            "btopt", "K", BcArrayFloatState, "one-dimensional",
            "`(xnum, 2)` : Temperature boundary condition along top boundary: \n"
            *"\t tk[1][j] = btopt[j][1] + tk[2][j]*btop[j][2]",
        ),
        bbottomt = ArrayData(
            "bbottomt", "K", BcArrayFloatState, "one-dimensional",
            "`(xnum, 2)` : Temperature boundary condition along bottom boundary: \n"
            *"\t tk[ynum][j] = bbottomt[j][1] + tk[ynum-1][j]*bbottomt[j][2]",
        ),
        bleftt = ArrayData(
            "bleftt", "K", BcArrayFloatState, "one-dimensional",
            "`(ynum, 2)` : Temperature boundary condition along left boundary: \n"
            *"\t tk[i][1] = bleftt[i][1] + bleftt[i][2]*tk[i][2]",
        ),
        brightt = ArrayData(
            "brightt", "K", BcArrayFloatState, "one-dimensional",
            "`(ynum, 2)` : Temperature boundary condition along right boundary: \n"
            *"\t tk[i][xnum] = brightt[i][1] + brightt[i][2]*tk[i][xnum-1]",
        ),

        # Component-wise Velocity Boundary Conditions Arrays
        btopx = ArrayData(
            "btopx", "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vx boundary condition along top boundary: \n"
            *"\t vx[1,j,:] = btopx[j,1] + vx[2,j,:]*btopx[j,2]",
        ),
        btopy = ArrayData(
            "btopy", "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vy boundary condition along top boundary: \n"
            *"\t vy[1,j,:] = btopy[j,1] + vy[2,j,:]*btopy[j,2]",
        ),
        btopz = ArrayData(
            "btopz", "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vz boundary condition along top boundary: \n"
            *"\t vz[1,j,:] = btopz[j,1] + vz[2,j,:]*btopz[j,2]",
        ),
        bbottomx = ArrayData(
            "bbottomx", "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vx boundary condition along bottom boundary: \n"
            *"\t vx[ynum+1,j,:] = bbottomx[j,1] + vx[ynum, j,:]*bbottomx[j,2]",
        ),
        bbottomy = ArrayData(
            "bbottomy", "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vy along bottom boundary: \n"
            *"\t vy[ynum,j,:] = bbottomy[j,1] + vy[ynum-1, j,:]*bbottomy[j,2]",
        ),
        bbottomz = ArrayData(
            "bbottomz", "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vz boundary condition along bottom boundary: \n"
            *"\t vz[ynum+1,j,:] = bbottomz[j,1] + vz[ynum, j,:]*bbottomz[j,2]",
        ),
        bleftx = ArrayData(
            "bleftx", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vx boundary condition along left boundary: \n"
            *"\t vx[i,1] = bleftx[i,1] + vx[i, 2]*bleftx[i,2]",
        ),
        blefty = ArrayData(
            "blefty", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vy boundary condition along left boundary: \n"
            *"\t vy[i,1] = blefty[i,1] + vy[i, 2]*blefty[i,2]",
        ),
        bleftz = ArrayData(
            "bleftz", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vz boundary condition along left boundary: \n"
            *"\t vz[i,1,:] = bleftz[i,1] + vz[i,2,:]*bleftz[i,2]",
        ),
        brightx = ArrayData(
            "brightx", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vx boundary condition along right boundary: \n"
            *"\t vx[i,xnum] = brightx[i, 1] + vx[i, xnum-1]*brightx[i,2]",
        ),
        brighty = ArrayData(
            "brighty", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vy boundary condition along right boundary: \n"
            *"\t vy[i,xnum+1] = brighty[i,1] + vx[i, xnum]*brighty[i,2]",
        ),
        brightz = ArrayData(
            "brightz", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vz boundary condition along right boundary: \n"
            *"\t vz[i,xnum+1,:] = brightz[i,1] + vz[i,xnum,:]*brightz[i,2]",
        ),
        bfrontx = ArrayData(
            "bfrontx", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vx boundary condition along front boundary: \n"
            *"\t vx[i,j,1] = bfrontx[i,j,1] + vx[i,j,2]*bfrontx[i,j,2]",
        ),
        bfronty = ArrayData(
            "bfronty", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vy boundary condition along front boundary: \n"
            *"\t vy[i,j,1] = bfronty[i,j,1] + vy[i,j,2]*bfronty[i,j,2]",
        ),
        bfrontz = ArrayData(
            "bfrontz", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vz boundary condition along front boundary: \n"
            *"\t vz[i,j,1] = bfrontz[i,j,1] + vz[i,j,2]*bfrontz[i,j,2]",
        ),
        bbackx = ArrayData(
            "bbackx", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vx boundary condition along back boundary: \n"
            *"\t vx[i,j,znum+1] = bbackx[i,j,1] + vx[i,j,znum]*bbackx[i,j,2]",
        ),
        bbacky = ArrayData(
            "bbacky", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vy boundary condition along back boundary: \n"
            *"\t vy[i,j,znum+1] = bbacky[i,j,1] + vy[i,j,znum]*bbacky[i,j,2]",
        ),
        bbackz = ArrayData(
            "bbackz", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vz boundary condition along back boundary: \n"
            *"\t vz[i,j,znum  ] = bbackz[i,j,1] + vz[i,j,znum-1]*bbackz[i,j,2]",
        ),

        # Combined Velocity Boundary Conditions Arrays
        btop = ArrayData(
            "btop", "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 6)` : Velocity boundary conditions along top boundary: \n"
            *"\t vx[1,j,:] = btop[j,1] + vx[2,j,:]*btop[j,2] \n"
            *"\t vy[1,j,:] = btop[j,3] + vy[2,j,:]*btop[j,4] \n"
            *"\t vz[1,j,:] = btop[j,5] + vz[2,j,:]*btop[j,6]"
        ),
        bbottom = ArrayData(
            "bbottom", "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 6)` : Velocity boundary conditions along bottom boundary: \n"
            *"\t vx[ynum+1,j,:] = bbottom[j,1] + vx[ynum  ,j,:]*bbottom[j,2] \n"
            *"\t vy[ynum  ,j,:] = bbottom[j,3] + vy[ynum-1,j,:]*bbottom[j,4] \n"
            *"\t vz[ynum+1,j,:] = bbottom[j,5] + vz[ynum  ,j,:]*bbottom[j,6]",
        ),
        bleft = ArrayData(
            "bleft", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 6)` : Velocity boundary conditions along left boundary: \n"
            *"\t vx[i,1,:] = bleft[i,1] + vx[i,2,:]*bleft[i,2] \n"
            *"\t vy[i,1,:] = bleft[i,3] + vy[i,2,:]*bleft[i,4] \n"
            *"\t vz[i,1,:] = bleft[i,5] + vz[i,2,:]*bleft[i,6]",
        ),
        bright = ArrayData(
            "bright", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 6)` : Velocity boundary conditions along right boundary: \n"
            *"\t vx[i,xnum  ,:] = bright[i,1] + vx[i,xnum-1,:]*bright[i,2] \n"
            *"\t vy[i,xnum+1,:] = bright[i,3] + vy[i,xnum  ,:]*bright[i,4] \n"
            *"\t vz[i,xnum+1,:] = bright[i,5] + vz[i,xnum  ,:]*bright[i,6]",
        ),
        bfront = ArrayData(
            "bfront", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 6)` : Velocity boundary conditions along front boundary: \n"
            *"\t vx[i,j,1] = bfront[i,1] + vx[i,j,2]*bfront[i,2] \n"
            *"\t vy[i,j,1] = bfront[i,3] + vy[i,j,2]*bfront[i,4] \n"
            *"\t vz[i,j,1] = bfront[i,5] + vz[i,j,2]*bfront[i,6]",
        ),
        bback = ArrayData(
            "bback", "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 6)` : Velocity boundary conditions along back boundary: \n"
            *"\t vx[i,j,znum+1] = bback[i,1] + vx[i,j,znum]*bback[i,2] \n"
            *"\t vy[i,j,znum+1] = bback[i,3] + vy[i,j,znum]*bback[i,4] \n"
            *"\t vz[i,j,znum  ] = bback[i,5] + vz[i,j,znum-1]*bback[i,6]",
        ),
    )
end

end # module