# Building and Executing Models with the EarthBox API

EarthBox models can be constructed and executed using an API. The recommended
way to do this is by calling API functions (see [EarthBox API Overview](@ref)) from a 
script called [`Model.jl`](@ref) that includes a module called `Model` containing the 
API calls and a function called `run_case()` that can be executed via command 
line. A script called [`Plot.jl`](@ref) is built by the user for plotting output 
files using customizable functions for scalar grids, markers, yield strength profiles
and convergence of the non-linear Picard iterations. 

The following scripts are commonly used to build, initialize and execute
EarthBox models using the API:
- [`Model.jl`](@ref)
- [`Materials.jl`](@ref)
- [`Plot.jl`](@ref)

!!! note "`run_earthbox()`"
    The [`Model.jl`](@ref) and [`Plot.jl`](@ref) scripts are required by the 
    `run_earthbox()` function used for sequential background command-line execution 
    (see [Background Command-line Execution](@ref)) and for execution on remote 
    machines (see [Batch Remote Execution](@ref)).

# `Model.jl`

The `Model.jl` script defines key paths, executes initialization API functions
and runs model time steps. The example below is for an extensional model using 
a weak fault zone to localize deformation. Material input is defined 
by importing a function from a script called [Materials.jl](@ref) that defines 
the a dictionary called `materials_input_dict` with types, domains, colors and 
material names from a material library where properties are defined. Alternatively,
material input could be defined using a yamal formatted file (See 
[MarkerMaterials](@ref) and [Material Input](@ref) for more details).

The `Model.jl` script should be setup to manage model case names using EarthBox 
API functions (see [EarthBox.GetArgs.get_case_name_from_cl_args]). The format 
used for model cases is `case0`, `case1`, `case2`, etc (default=`case0`). Model 
output is sent to the model output directory named `<case_name>_output` where 
`<case_name>` is the name of the model case. 

The model output path is defined using the `case_name` and the `ROOT_PATH_OUTPUT`
constant unless the `model_output_path` is provided via command line argument.

The `Model.jl` can be used with the module `CaseInputs.jl` to define model input 
parameters based on the model case name (see [Case Management](@ref)). This is 
accomplished by including the `CaseInputs.jl` module and calling the 
`define_case_inputs` function with the `case_name` as an argument.

## Command-line Arguments for `Model.jl`
- `case_name=<case_name>`
    - The name of the model case like "case0", "case1", "case2", etc. Default is 
        "case0". The case name will be used to define inputs if this script is
        integrated with the `CaseInputs.jl` module. The case name is also used to
        define the model output path if the model output path is not provided via a 
        command line argument. If the case name is passed as an argument to the 
        `run_case` function, the `case_name` command line argument will be ignored.

-`model_output_path=<model_output_path>`
    - The path to the model output directory. If not provided, then the model output
       path is defined using the `case_name` and the `ROOT_PATH_OUTPUT`.

## Using `Model.jl` From the REPL

```julia
include("Model.jl")
Model.run_case(case_name=<case_name>)
```

If no `<case_name>` is provided, then the default case name is `case0`.

## Using `Model.jl` From Command-line

```bash
julia Model.jl case_name=<case_name>
```

where `<case_name>` is the name of the model case and the model output path is 
defined automatically using the `case_name` and the constant `ROOT_PATH_OUTPUT`.

You can also specify the model output path via a command line argument:

```bash
julia Model.jl case_name=<case_name> model_output_path=<model_output_path>
```

## Example `Model.jl`: Extension with Weak Fault

```julia
module Model

using EarthBox
include("Materials.jl")
import .Materials: get_materials_input_dict

const ROOT_PATH_OUTPUT = "/mnt/extradrive1/earthbox_output/extension_to_sfs/extension_strong_zones"
const MATERIAL_COLLECTION_PATH = MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1.path

# Get the input parameters object so names can be accessed programmatically
const PARAMS = get_eb_parameters()

# Constants are defined for parameters used in multiple initialization functions
const adiabatic_gradient       = 0.4 # K/km
const xsize                    = 500_000.0 # meters
const ysize                    = 160_000.0 # meters
const thick_air                = 10_000.0 # meters
const thick_crust              = 32_000.0 # meters
const thick_upper_crust        = 22_000.0 # meters
const thick_lith               = 125_000.0 # meters
const marker_spacing           = 100.0 # meters (50 m for high resolution case)
const grid_spacing_high_res    = 500.0 # meters (200 m for high resolution case)
const avg_grid_spacing_low_res = 2000.0 # meters
const temperature_base_lith_celsius = 1330.0

function run_case(;case_name::String="case0")::Nothing
    print_info("Running model for case: $case_name")
    # model output directories are named <case_name>_output where <case_name> is 
    # the name of the model case
    model_output_path = get_model_output_path(case_name, ROOT_PATH_OUTPUT)
    eb = setup_model(model_output_path=model_output_path)
    PRINT_SETTINGS.print_performance = false
    # Here we specify model duration and let EarthBox calculate the viscoelastic
    # time step based on minimum grid spacing and extension rate. The number of
    # time steps is estimated based on the `model_duration_myr` parameter, velocity 
    # stepping and the calculated viscoelastic time step. Alternatively, you could specify 
    # `timestep_viscoelastic` in seconds and `ntimestep_max`. This will deactivate
    # the automated calculations and requires a bit of math to define adequate values.
    # If the calculated velocity is higher than the half-spreading rate leading to
    # a reduced time step , you may may have to set the model duration to a
    # larger value to achieve the desired actual model duration.
    run_time_steps(
        eb,
        make_backup           = true,
        model_duration_myr    = 50.0,
        timestep_out          = ConversionFuncs.years_to_seconds(500_000.0)
    )
    return nothing
end

function setup_model(;model_output_path::String)::EarthBoxState
    # Number of basic nodes in x- and y-directions (xnum, ynum) are calculated 
    # from the T-type refinement parameters and used to initialize the staggered 
    # grid. Alternatively, you could define xnum and ynum directly as input 
    # parameters and provided t-type parameters during Staggered grid initialization.
    # However, this will require putting some thought into defining a sufficient
    # number of basic nodes in x- and y-directions.
    ttype_refinement_parameters = Dict(
        PARAMS.xo_highres.name => 150_000.0, # meters
        PARAMS.xf_highres.name => 350_000.0, # meters
        PARAMS.yf_highres.name => 50_000.0, # meters
        PARAMS.dx_highres.name => grid_spacing_high_res, # meters (200 m for high resolution case)
        PARAMS.dx_lowres.name  => avg_grid_spacing_low_res, # meters
        PARAMS.dy_highres.name => grid_spacing_high_res, # meters (200 m for high resolution case)
        PARAMS.dy_lowres.name  => avg_grid_spacing_low_res # meters
    )
    eb = EarthBoxState(
        restart_from_backup         = false,
        xsize                       = xsize, # meters
        ysize                       = ysize, # meters
        dx_marker                   = marker_spacing, # meters (50 m for high resolution case)
        dy_marker                   = marker_spacing, # meters (50 m for high resolution case)
        ttype_refinement_parameters = ttype_refinement_parameters,
        use_mumps                   = false,
        nprocs                      = 8,
        use_internal_mumps          = true,
        paths                       = Dict("output_dir" => model_output_path)
    )
    initialize_model_input(eb)
    return eb
end

function initialize_model_input(eb::EarthBoxState)::Nothing
    initialize_staggered_grid(eb)
    initialize_geometry(eb)
    initialize_boundary_conditions(eb)
    initialize_marker_coordinates(eb)
    initialize_marker_materials(eb)
    initialize_marker_temperature(eb)
    initialize_marker_plasticity(eb)
    initialize_rock_properties(eb)
    initialize_stokes_continuity_solver(eb)
    initialize_heat_solver(eb)
    initialize_advection_model(eb)
    initialize_interpolation_model(eb)
    initialize_surface_processes_model(eb)
    return nothing
end

function initialize_staggered_grid(eb::EarthBoxState)::Nothing
    # T-type refinement parameters are used from EarthBoxState initialization
    # and not required to be set here. However, if T-type refinement parameters
    # were not previously defined, they must be set here.
    StaggeredGrid.initialize!(
        eb.model_manager.model, grid_type = grid_type_names.TtypeRefinedGrid)
    return nothing
end

function initialize_geometry(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    MaterialGeometry.StickyAirGeometry.initialize!(model, thick_air=thick_air)
    MaterialGeometry.EarthLayering.initialize!(
        model,
        thick_lith        = thick_lith, # meters
        thick_crust       = thick_crust, # meters
        thick_upper_crust = thick_upper_crust, # meters
    )
    # The strong zone boundaries are set using T-type refinement parameters
    # and not required to be set here. However, if T-type refinement parameters
    # were not previously defined, `x_left_strong` and `x_right_strong` must be set 
    # here.
    MaterialGeometry.LithoStrongZones.initialize!(model, iuse_strong_zones = 1)
    return nothing
end

function initialize_boundary_conditions(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    BoundaryConditions.initialize!(
        model, model_type = bc_model_type_names.LithosphericExtensionFixedBoundaries)
    BoundaryConditions.Pressure.initialize!(model, pressure_bc=1e5)
    BoundaryConditions.Temperature.initialize!(
        model,
        temperature_top    = ConversionFuncs.celsius_to_kelvin(0.0),
        temperature_bottom = ConversionFuncs.celsius_to_kelvin(
            temperature_base_lith_celsius 
            + adiabatic_gradient*(ysize - thick_air - thick_lith)/1000.0
        )
    )
    BoundaryConditions.Velocity.initialize!(
        model,
        full_velocity_extension = ConversionFuncs.cm_yr_to_m_s(0.1),
        full_velocity_extension_step1 = ConversionFuncs.cm_yr_to_m_s(0.5),
        )
    # It is often required to start extension models with a relatively low extension rate
    # and then later increase the extension rate to avoid numerical issues that
    # may manifest themselves as erroneous high strain zones. This issue is
    # overcome by using velocity stepping. In this example we turn on velocity stepping,
    # specify and step time and let EarthBox calculate the the velocity step factor 
    # (`velocity_step_factor`) and time step adjustment factor (`timestep_adjustment_factor`).
    # based on the extension velocities specified above (`full_velocity_extension` and 
    # `full_velocity_extension_step1`). Alternatively, you could specify the velocity step
    # factor (`velocity_step_factor`) and time step adjustment factor (`timestep_adjustment_factor`) 
    # directly.
    BoundaryConditions.VelocityStep.initialize!(
        model,
        iuse_velocity_step = 1,
        velocity_step_time = 0.5, # Myr
    )
    return nothing
end
 
function initialize_marker_coordinates(eb::EarthBoxState)::Nothing
    Markers.MarkerCoordinates.initialize!(
        eb.model_manager.model, 
        marker_distribution=marker_distribution_names.Randomized
        )
    return nothing
end

function initialize_marker_materials(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    Markers.MarkerMaterials.initialize!(
        model,
        material_model       = material_model_names.LithosphericExtensionLateralStrongZones,
        paths                = Dict("materials_library_file" => MATERIAL_COLLECTION_PATH),
        materials_input_dict = get_materials_input_dict(),
        viscosity_min        = 1e18, # Pa.s
        viscosity_max        = 1e26 # Pa.s
    )
    Markers.MarkerMaterials.MarkerStressLimits.initialize!(
        model, 
        powerlaw_stress_min = 1e4, # Pa
        yield_stress_min    = 0.0, # Pa
        yield_stress_max    = 1e32 # Pa
    )
    return nothing
end

function initialize_marker_temperature(eb::EarthBoxState)::Nothing
    Markers.MarkerTemperature.initialize!(
        eb.model_manager.model,
        initial_temperature_model=initial_temperature_names.AnalyticalThreeLayer,
        parameters=Dict(
            PARAMS.temperature_base_lith.name       => (
                ConversionFuncs.celsius_to_kelvin(temperature_base_lith_celsius)),
            PARAMS.amplitude_perturbation.name      => 4000.0, # meters
            PARAMS.width_perturbation.name          => 10_000.0, # meters
            PARAMS.thick_thermal_lithosphere.name   => thick_lith, # meters
            PARAMS.adiabatic_gradient.name          => adiabatic_gradient, # K/km
            PARAMS.conductivity_upper_crust.name    => 2.25, # W/m/K
            PARAMS.conductivity_lower_crust.name    => 2.2, # W/m/K
            PARAMS.conductivity_mantle.name         => 2.0, # W/m/K
            PARAMS.heat_production_upper_crust.name => 1.25e-6, # W/m^3
            PARAMS.heat_production_lower_crust.name => 0.2e-6, # W/m^3
            PARAMS.heat_production_mantle.name      => 0.0, # W/m^3
        )
    )
    return nothing
end

function initialize_marker_plasticity(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    Markers.MarkerFriction.initialize!(
        model, initialization_model = friction_init_names.Regular)
    Markers.MarkerCohesion.initialize!(model)
    Markers.MarkerPreexponential.initialize!(model)
    return nothing
end

function initialize_rock_properties(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    RockProperties.initialize!(
        model,
        density_model              = density_model_names.Liao14,
        thermal_conductivity_model = thermal_conductivity_model_names.SekiguchiWaples,
        rhocp_model                = rhocp_model_names.TemperatureDependentWaples,
        iuse_sed_porosity          = 1
    )
    return nothing
end

function initialize_stokes_continuity_solver(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    StokesContinuitySolver.initialize!(
        model,
        velocity_type                = velocity_type_names.VelocityFromStokesSolver,
        gravity_y                    = 9.8, # m/s^2
        iuse_interface_stabilization = 1
    )
    GlobalPlasticityLoop.initialize!(
        model,
        global_plasticity_loop = global_plasticity_names.NodalPlasticityLoop,
        tolerance_picard       = 1e-2,
        nglobal                = 3
    )
    return nothing
end

function initialize_heat_solver(eb::EarthBoxState)::Nothing
    HeatSolver.initialize!(
        eb.model_manager.model,
        iuse_heat              = 1,
        iuse_adiabatic_heating = 1,
        iuse_shear_heating     = 1,
        iuse_sticky_correction = 1,
        # The max_temp_change is used for the adaptive time stepping scheme used
        # to solve the heat equation.
        max_temp_change        = 70.0 # K
    )
    return nothing
end

function initialize_advection_model(eb::EarthBoxState)::Nothing
    Advection.initialize!(
        eb.model_manager.model,
        advection_scheme         = advection_scheme_names.RungeKutta4thOrder,
        # Decreasing the marker_cell_displ_max can increase numerical accuracy.
        marker_cell_displ_max    = 0.75, # fraction
        subgrid_diff_coef_temp   = 1.0,
        subgrid_diff_coef_stress = 1.0
    )    
    return nothing
end

function initialize_interpolation_model(eb::EarthBoxState)::Nothing
    Interpolation.initialize!(
        eb.model_manager.model,
        iuse_harmonic_avg_normal_viscosity = 1,
        iuse_initial_temp_interp           = 1
    )
    return nothing
end

function initialize_surface_processes_model(eb::EarthBoxState)::Nothing
    model = eb.model_manager.model
    
    SurfaceProcesses.Topography.initialize!(
        model,
        iuse_topo            = 1,
        node_advection       = topo_node_advection_names.RungeKuttaWithInterp,
        dx_topo              = 200.0, # meters
        topo_xsize           = xsize, # meters
        nsmooth_top_bottom   = 2,
        marker_search_factor = 2.0    
    )
    
    SurfaceProcesses.Sealevel.initialize!(
        model,
        option_name               = sealevel_option_names.AveragePressure,
        y_sealevel                = thick_air, # meters
        base_level_shift_end_time = 16.0, # Myr
        base_level_shift          = 0.0, # meters
    )

    SurfaceProcesses.Sealevel.RelativeBaseLevel.initialize!(
        model,
        iuse_linear_segments                   = 0,
        thickness_upper_continental_crust_ref = thick_upper_crust, # meters
        thickness_lower_continental_crust_ref = thick_crust - thick_upper_crust, # meters
        thickness_lithosphere_ref             = thick_lith, # meters
        gridy_spacing_ref                     = 100.0, # meters
        temperature_top_ref                   = 0.0, # Celsius
        temperature_moho_ref                  = 600.0, # Celsius
        temperature_base_lith_ref             = temperature_base_lith_celsius, # Celsius
        adiabatic_gradient_ref                = adiabatic_gradient, # K/km
    )
    return nothing

end

function main()::Nothing
    run_case(case_name=get_case_name_from_cl_args())
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Model.main()
end

```

# `Materials.jl`

The `Materials.jl` script is used to define the material inputs used in the model.
This script is not required as material input dictionaries could be built in
`Model.jl` or material input could be defined using a yamal file (See 
[MarkerMaterials](@ref) and [Material Input](@ref) for more details).

```julia
module Materials

using EarthBox

function get_materials_input_dict()::MaterialsDictType
    mat_names = MaterialLibrary().lithospheric_deformation.lithospheric_deformation_eb1.materials
    types = MaterialTypesRegistry()
    domains = MaterialDomainsRegistry()

    return MaterialsDictType(
        # Sticky Domain
        Int16(1) => MaterialDictType(
            "mat_name" => mat_names.sticky_air,
            "mat_type" => types.sticky_air,
            "mat_domain" => domains.atmosphere,
            "red_fraction" => 255/255,
            "green_fraction" => 255/255,
            "blue_fraction" => 255/255,
        ),
        Int16(2) => MaterialDictType(
            "mat_name" => mat_names.sticky_water,
            "mat_type" => types.sticky_water,
            "mat_domain" => domains.ocean,
            "red_fraction" => 0/255,
            "green_fraction" => 255/255,
            "blue_fraction" => 255/255,
        ),
        # Sedimentary Basin
        Int16(3) => MaterialDictType(
            "mat_name" => mat_names.clastic_sediment,
            "mat_type" => types.sediment,
            "mat_domain" => domains.sedimentary_basin,
            "red_fraction" => 0.89411765,
            "green_fraction" => 0.58039216,
            "blue_fraction" => 0.28627451,
        ),
        # Felsic Continental Crust
        Int16(4) => MaterialDictType(
            "mat_name" => mat_names.felsic_continental_crust,
            "mat_type" => types.general,
            "mat_domain" => domains.upper_continental_crust,
            "red_fraction" => 255/255,
            "green_fraction" => 153/255,
            "blue_fraction" => 153/255,
        ),
        Int16(5) => MaterialDictType(
            "mat_name" => mat_names.felsic_continental_crust_strong_zone,
            "mat_type" => types.general,
            "mat_domain" => domains.upper_continental_crust_strong_zone,
            "red_fraction" => 255/255,
            "green_fraction" => 163/255,
            "blue_fraction" => 163/255,
        ),
        # Mafic Continental Crust
        Int16(6) => MaterialDictType(
            "mat_name" => mat_names.mafic_continental_crust,
            "mat_type" => types.general,
            "mat_domain" => domains.lower_continental_crust,
            "red_fraction" => 255/255,
            "green_fraction" => 200/255,
            "blue_fraction" => 200/255,
        ),
        Int16(7) => MaterialDictType(
            "mat_name" => mat_names.mafic_continental_crust_strong_zone,
            "mat_type" => types.general,
            "mat_domain" => domains.lower_continental_crust_strong_zone,
            "red_fraction" => 255/255,
            "green_fraction" => 210/255,
            "blue_fraction" => 210/255,
        ),
        # Continental Mantle Lithosphere
        Int16(8) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_continental_lithosphere,
            "mat_type" => types.general,
            "mat_domain" => domains.mantle_lithosphere,
            "red_fraction" => 0.0,
            "green_fraction" => 153/255,
            "blue_fraction" => 153/255,
        ),
        Int16(9) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_continental_lithosphere_strong_zone,
            "mat_type" => types.general,
            "mat_domain" => domains.lithospheric_mantle_strong_zone,
            "red_fraction" => 0.0,
            "green_fraction" => 163/255,
            "blue_fraction" => 163/255,
        ),
        # Asthenosphere
        Int16(10) => MaterialDictType(
            "mat_name" => mat_names.ultramafic_asthenosphere_dry_fertile,
            "mat_type" => types.general,
            "mat_domain" => domains.asthenosphere,
            "red_fraction" => 0.0,
            "green_fraction" => 200/255,
            "blue_fraction" => 153/255,
        )
    )
end

end # module
```


# `Plot.jl`

The `Plot.jl` script is a customized plotting module with user defined plotting 
functions and command-line plotter that plots output files stored in the model 
output directory. Note that model output directories are named `<case_name>_output` 
where `<case_name>` is the name of the model case like `case0`, `case1`, `case2`, etc (default=`case0`).

For details on command-line arguments and examples, see [Command-line Plotter](@ref).

## Using `Plot.jl` From Command-line:

```bash
julia Plot.jl <plot_option_name> case_name=<case_name> istart=1 iend=100
```

where `istart` and `iend` are the starting and ending time steps, `<case_name>` 
is the name of the model case and `<plot_option_name>` is the name of the plot 
option that can be set equal to one of the following values:

- `marker_plots`: Plot markers.
- `scalar_plots`: Plot scalars.
- `velocity_plots`: Plot velocity.
- `stokes_convergence_plots`: Plot Stokes convergence.
- `yield_strength_plots`: Plot yield strength.

## Using `Plot.jl` From REPL

```julia
include("Plot.jl")
Plot.marker_plots(
    model_output_path = "/path/to/model_output",
    istart = 1,
    iend = 100
);
```

!!! note "istart and iend"
    `istart` and `iend` are not used by the stokes convergence and yield 
    strength plotting functions.

## Example `Plot.jl`: Extension with Weak Fault

```julia
module Plot

using EarthBox

include("Model.jl")
include("Materials.jl")
import .Model: MATERIAL_COLLECTION_PATH, ROOT_PATH_OUTPUT
import .Materials: get_materials_input_dict 

const dimensions = (0.0, 500.0, 0.0, 160.0)
const xyspacing = (50.0, 10.0)
const model_figsize = (10.0, 5.0)

function scalar_plots(;
    model_output_path::String,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    plot_scalars(
        scalar_name = :temperature,
        model_output_path = model_output_path,
        dimensions = dimensions,
        xyspacing = xyspacing,
        istart = istart,
        iend = iend,
        figsize = model_figsize
    )
    plot_scalars(
        scalar_name = :viscosity,
        model_output_path = model_output_path,
        dimensions = dimensions,
        xyspacing = xyspacing,
        istart = istart,
        iend = iend,
        figsize = model_figsize
    )
    return nothing
end

function velocity_plots(;
    model_output_path::String,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    plot_velocity(
        model_output_path = model_output_path,
        dimensions        = dimensions,
        xyspacing         = xyspacing,
        istart            = istart,
        iend              = iend,
        decimation_factor = 8,
        scale_factor      = 10.0,
        figsize = model_figsize
    )
    return nothing
end

function marker_plots(;
    model_output_path::String,
    istart::Int64 = 1,
    iend::Union{Int64, Nothing} = nothing
)::Nothing
    plot_markers(
        plot_type=:CompositionHeatFlowGravity,
        model_output_path=model_output_path,
        material_library_file_path=MATERIAL_COLLECTION_PATH,
        materials_input_dict=get_materials_input_dict(),
        # General plotting parameters
        dimensions=dimensions,
        xyspacing=xyspacing,
        istart=istart, iend=iend,
        figsize = model_figsize,
        xy_location_contour_legend = (1.0, 5.0),
        text_box_font_size = 5,
        # General marker parameters
        # (reducing the decimation factor reduces pixelation of the plot)
        marker_size=2.0, decimation_factor=5,
        plot_mesh=0, mesh_line_width=0.1,
        plot_contour_labels=1, contour_line_width=1.0, contour_line_color="black",
        # Topography parameters
        plot_topography=1, topo_line_width=1.0, 
        topo_line_color="red",
        # Base level parameters
        plot_base_level=1, base_level_line_width=1.0,
        base_level_line_color="blue",
        # Plastic strain parameters
        plot_plastic_strain = true,
        strain_min = 1.0, strain_max = 6.0,
        strain_contour_interval = 0.25, strain_cmap="inferno",
        # Plastic strain rate parameters
        plot_plastic_strain_rate = 1,
        strain_rate_min = -18, strain_rate_max = -12,
        strain_rate_contour_interval = 0.5, strain_rate_cmap = "Reds",
        # Temperature contours parameters
        plot_temperature_contours=1,
        temperature_min=100.0, temperature_max=1400.0,
        temperature_contour_interval=100.0,
        temperature_number_format="%6.0f",
        temperature_contour_color="black",
        # Heat flow parameters
        heatflow_min=0.0, heatflow_max=200.0,
        heatflow_spacing=50.0,
        # Gravity Parameters (mgal)
        gravity_min=-1000.0, gravity_max=1000.0,
        gravity_spacing=200.0
    )
    return nothing
end

function yield_strength_plots(;
    model_output_path::String
)::Nothing
    plot_yield_strength(
        model_output_path=model_output_path,
        material_library_file_path=MATERIAL_COLLECTION_PATH,
        materials_input_dict=get_materials_input_dict(),
        depth_plot_limit=150.0,
        depth_axis_spacing=10.0,
        temperature_plot_limit=1400.0,
        temperature_axis_spacing=200.0,
        yield_stress_plot_limit=600.0,
        yield_stress_axis_spacing=50.0,
        strain_rate=1e-15,
        thickness_upr_cont_crust_meters=22_000.0,
        thickness_lwr_cont_crust_meters=10_000.0,
        thickness_lithosphere_meters=125_000.0,
        thickness_asthenosphere_meters=25_000.0,
        dy_meters=100.0,
        expansivity=3e-5,
        compressibility=1e-11,
        density_upper_continental_crust=2822.72,
        density_lower_continental_crust=2900.0,
        density_mantle_lithosphere=3250.0,
        density_asthenosphere=3300.0,
        iuse_linear_segments=false,
        temperature_top_celsius=0.0,
        temperature_moho_celsius=600.0,
        temperature_base_lith_celsius=1330.0,
        adiabatic_gradient_kelvin_km=0.4,
        conductivity_upper_crust=2.25,
        conductivity_lower_crust=2.2,
        conductivity_mantle=2.0,
        heat_production_upper_crust=1.8e-6,
        heat_production_lower_crust=0.5e-6,
        heat_production_mantle=0.0,
        thickness_thermal_lithosphere=125_000.0,
        iuse_fluid_pressure_for_yield=1,
        log10_viscosity_plot_min=20.0,
        log10_viscosity_plot_max=26.0,
        extension=".pdf",
        figsize=(3.0, 6.0)
        )
    return nothing
end

function stokes_convergence_plots(;
    model_output_path::String
)::Nothing
    plot_stokes_convergence(
        model_output_path=model_output_path,
        figsize=(16, 4),
        xmin=0,
        xmax=200,
        xspacing=50,
        log_l2_ymin=-7,
        log_l2_ymax=1,
        log_l2_yspacing=1,
        plot_title="TestCase",
        axis_label_size=8,
        axis_title_size=8,
        legend_font_size=8,
        xtick_label_size=6,
        ytick_label_size=6,
        annotation_font_size=4
        )
    return nothing
end

function main()::Nothing
    run_cl_plotter(
        root_path_output              = ROOT_PATH_OUTPUT,
        marker_plots_func             = marker_plots,
        scalar_plots_func             = scalar_plots,
        velocity_plots_func           = velocity_plots,
        stokes_convergence_plots_func = stokes_convergence_plots,
        yield_strength_plots_func     = yield_strength_plots
    )
    return nothing
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    Plot.main()
end

```

