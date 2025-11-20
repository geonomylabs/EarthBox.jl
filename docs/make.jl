using Documenter, EarthBox
using DocumenterCitations

path = joinpath(@__DIR__, "src", "bibliography", "bibliography.bib")
bib = CitationBibliography(path, style=:authoryear)

introduction = "Introduction" => "index.md"

install_pages = [
    "Installation" => "installation/installation.md",
    "Using System MPI Binaries" => "installation/MPI_system.md",
    "Using Custom MUMPS Libraries" => "installation/MUMPS_custom.md"
    ]

model_pages = [
    "Model Management with the EarthBox API" => [
        "Model Building, Execution and Plotting" => "models/model_execution/api_building_and_execution.md",
        "Using Input Files" => "models/model_execution/input_file_execution.md",
        "Background Command-line Execution" => "models/model_execution/background_command_line_execution.md",
        "Batch Remote Execution" => "models/model_execution/remote_execution.md",
        "Batch Local Execution" => "models/model_execution/batch_local_execution.md",
        "Plotting Marker Swarms in Parallel" => "models/model_execution/parallel_marker_plotting.md"
    ],
    "Models" => [
        "Model List" => "models/model_list.md",
        "Polyphase Extension n17" => "models/lithospheric_extension/polyphase_extension_n17.md",
        "Extension with Weak Fault" => "models/lithospheric_extension/extension_weak_fault.md",
        "Extension with Softening b14" => "models/lithospheric_extension/extension_with_softening_b14.md",
        "SDR Formation" => "models/sdr_formation/sdr_formation.md",
        "Visco-elasto-plastic Extension with Downhill Diffusion" => "models/vep_extension_downhill_diffusion/vep_extension_downhill_diffusion.md",
        "Slab Retreat" => "models/subduction/slab_retreat.md",
        "Sandbox Extension" => "models/sandbox/sandbox_extension.md",
        "Sandbox Shortening" => "models/sandbox/sandbox_shortening.md",
        "Spreading Center" => "models/spreading_center/spreading_center.md",
        "Hot Box Melting" => "models/hot_box_melting/hot_box_melting.md",
        "Magmatic Crust Formation: Central Atlantic" => "models/regional/catl.md",
        "Magma Poor Margin: Iberia" => "models/regional/iberia.md",
    ]
]

methods_pages = [
    "Introduction" => "methods/introduction.md",
    "Governing Equations" => [
        "Conservation of Momentum and Mass" => "methods/governing_equations/momentum_and_mass.md",
        "Conservation of Energy" => "methods/governing_equations/energy.md"
    ],
    "Discretization" => [
        "Staggered Grid Structure" => "methods/discretization/grid_structure.md",
        "Staggered Grid Generation" => "methods/discretization/grid_generation.md",
        "Marker Distribution" => "methods/discretization/marker_distribution.md",
        "Discretization of the 2D Stokes-continuity Equations" => "methods/discretization/stokes2D.md",
        "Discretization of Heat-conduction Equation" => "methods/discretization/conduction.md",
    ],
    "Boundary Conditions" => [
        "Stokes Boundary Conditions" => "methods/discretization/stokes_bcs.md",
        "Heat Boundary Conditions" => "methods/discretization/heat_bcs.md"
    ],
    "Marker-Grid Interpolation" => "methods/interpolation/marker_grid_interpolation.md",
    "Stokes-Continuity Algorithm" => "methods/stokes_algorithm/stokes_algorithm.md",
    "Heat Algorithm" => "methods/heat_algorithm/heat_algorithm.md",
    "Sediment Transport" => "methods/sediment_transport/sediment_transport.md",
    "Melt" => "methods/melt/melt.md",
    "Lava Flow" => "methods/lava_flow/lava_flow.md",
    "Viscous Creep" => "methods/viscous_creep/viscous_creep.md",
    "Strain Weakening" => "methods/strain_weakening/strain_weakening.md",
    "Hydrothermal Circulation" => "methods/hydrothermal/hydrothermal.md",
    "Serpentinization" => "methods/serpentinization/serpentinization.md"
]

thermo_mechanical_algorithm_pages = [
    "Main Steps" => "thermo_mechanical_algorithm/main_steps.md",
    "Pre-solver Steps" => [
        "thermo_mechanical_algorithm/pre_solver_steps/pre_solver_steps.md",
        "thermo_mechanical_algorithm/pre_solver_steps/sections/rock_property_update.md",
        "thermo_mechanical_algorithm/pre_solver_steps/sections/latent_heat_update.md",
        "thermo_mechanical_algorithm/pre_solver_steps/sections/hydrothermal_update.md",
        "thermo_mechanical_algorithm/pre_solver_steps/sections/marker_failure_prop_update.md",
        "thermo_mechanical_algorithm/pre_solver_steps/sections/marker_viscosity_update.md",
        "thermo_mechanical_algorithm/pre_solver_steps/sections/melt_damage_update.md",
    ],
    "Stokes-solver Steps" => "thermo_mechanical_algorithm/stokes_solver_steps/stokes_solver_steps.md",
    "Heat-solver Steps" => "thermo_mechanical_algorithm/heat_solver_steps/heat_solver_steps.md",
    "Advection-solver Steps" => "thermo_mechanical_algorithm/advection_solver_steps/advection_solver_steps.md",
    "Post-solver Steps" => [
        "Post-solver Steps" => "thermo_mechanical_algorithm/post_solver_steps/post_solver_steps.md",
        "Topography" => "thermo_mechanical_algorithm/post_solver_steps/sections/topography.md",
        "Melt Extraction" => "thermo_mechanical_algorithm/post_solver_steps/sections/melt_extraction.md",
        "Sediment Transport" => "thermo_mechanical_algorithm/post_solver_steps/sections/sediment_transport.md",
        "Lava Flow" => "thermo_mechanical_algorithm/post_solver_steps/sections/lava_flow.md",
        "Marker Transformation" => "thermo_mechanical_algorithm/post_solver_steps/sections/marker_transformation.md"
        ]
]

earthbox_api_pages = [
    "API Overview" => "api/earthbox_api.md",
    "EarthBoxState" => "api/earthbox_state/earthbox_state.md",
    "Model Management" => [
        "ModelManager" => "api/model_manager/model_manager.md",
        "XdmfTimeStepsManager" => "api/model_manager/xdmf/xdmf_timestep_manager.md"
        ],
    "Staggered Grid" => [
        "Initialization" => "api/staggered_grid/initialization.md"
    ],
    "Material Geometry" => [
        "StickyAirGeometry" => "api/material_geometry/sticky_air_geometry.md",
        "EarthLayering" => "api/material_geometry/earth_layering.md",
        "LithoStrongZones" => "api/material_geometry/litho_strong_zones.md",
        "Plume" => "api/material_geometry/plume.md",
        "CrustalHole" => "api/material_geometry/crustal_hole.md",
        "RayleighTaylor" => "api/material_geometry/rayleigh_taylor.md",
        "MobileWall" => "api/material_geometry/mobile_wall.md",
        "Sandbox" => "api/material_geometry/sandbox.md",
        "InternalVelocityZone" => "api/material_geometry/internal_velocity_zone.md",
        "FractureZone" => "api/material_geometry/fracture_zone.md",
        "WeakFault" => "api/material_geometry/weak_fault.md",
        "WeakSeed" => "api/material_geometry/weak_seed.md",
        ],
    "Boundary Conditions" => [
        "Initialization" => "api/boundary_conditions/boundary_conditions.md",
        "Pressure" => "api/boundary_conditions/bc_managers/pressure.md",
        "Temperature" => "api/boundary_conditions/bc_managers/temperature.md",
        "TransientBottomTemperature" => "api/boundary_conditions/bc_managers/transient_bottom_temperature.md",
        "Velocity" => "api/boundary_conditions/bc_managers/velocity.md",
        "VelocityFromStrainRate" => "api/boundary_conditions/bc_managers/velocity_from_strain_rate.md",
        "VelocityStop" => "api/boundary_conditions/bc_managers/velocity_stop.md",
        "VelocityStep" => "api/boundary_conditions/bc_managers/velocity_step.md",
        ],
    "Marker Initialization" => [
        "MarkerCoordinates" => "api/markers/coordinates.md",
        "MarkerTemperature" => "api/markers/temperature.md",
        "MarkerFriction" => "api/markers/friction.md",
        "MarkerCohesion" => "api/markers/cohesion.md",
        "MarkerDilatancy" => "api/markers/dilatancy.md",
        "MarkerViscousStrainSoftening" => "api/markers/viscous_strain_softening.md",
        "MarkerPreexponential" => "api/markers/preexponential.md",
        "MarkerStressLimits" => "api/markers/stress_limits.md",
        "MarkerBoundaryFriction" => "api/markers/boundary_friction.md",
    ],
    "Marker Materials" => [
        "MarkerMaterials" => "api/marker_materials/materials.md",
        "Material Collection Files" => "api/marker_materials/material_collection_files.md",
        "Material Properties" => "api/marker_materials/material_properties.md",
        "Material Input" => "api/marker_materials/material_input.md",
        "Material Types" => "api/marker_materials/material_types.md",
        "Material Domains" => "api/marker_materials/material_domains.md",
    ],
    "Rock Properties" => [
        "RockProperties" => "api/rock_properties/rock_properties.md",
        "Hydrothermal Circulation" => "api/rock_properties/hydrothermal/hydrothermal.md",
        "Serpentinization" => "api/rock_properties/serpentinization/serpentinization.md",
    ],
    "Stokes Continuity Solver" => [
        "StokesContinuitySolver" => "api/stokes_continuity_solver/stokes_continuity_solver.md",
        "GlobalPlasticityLoop" => "api/stokes_continuity_solver/global_plasticity_loop.md",
    ],
    "Heat Solver" => [
        "HeatSolver" => "api/heat_solver/heat_solver.md",
    ],
    "Marker Advection" => [
        "Advection" => "api/advection/advection.md",
    ],
    "Interpolation" => [
        "Interpolation" => "api/interpolation/interpolation.md",
    ],
    "Melt Model" => [
        "MeltModel" => "api/melt_model/melt_model.md",
        "MeltExtraction" => "api/melt_model/melt_extraction.md",
        "MeltExtrusion" => "api/melt_model/melt_extrusion.md",
        "MeltDamage" => "api/melt_model/melt_damage.md"
    ],
    "Surface Processes" => [
        "Topography" => "api/surface_processes/topography.md",
        "Sealevel" => "api/surface_processes/sealevel.md",
        "RelativeBaseLevel" => "api/surface_processes/relative_base_level.md",
        "SedimentTransport" => "api/surface_processes/sediment_transport.md",
        "SaltDeposition" => "api/surface_processes/salt.md",
    ],
    "Time Step Execution and Initialization" => [
        "Execution" => "api/run_time_steps/run_time_steps.md",
        "Initialization" => "api/run_time_steps/time_Step_initialization.md",
        "Examples" => "api/run_time_steps/examples.md",
    ],
    "Command-line Arguments" => "api/command_line_args/command_line_args.md",
    "Case Management" => "api/case_management/case_management.md",
]

material_library_pages = [
    "MaterialLibrary" => "material_library/material_library.md",
    "Collection: Lithospheric Deformation eb1" => "material_library/collections/lithospheric_deformation_eb1.md",
    "Collection: Benchmarks" => "material_library/collections/benchmarks.md",
    "Collection: Benchmarks SFS" => "material_library/collections/benchmarks_sfs.md",
    "Collection: Sandbox" => "material_library/collections/sandbox.md"
    ]

    plot_tools_pages = [
    "Command-line Plotter" => "plot_tools/cl_plotter.md",
    "Scalar Grid Plots" => "plot_tools/scalar_grid_plots.md",
    "Velocity Vector Plots" => "plot_tools/vector_plots.md",
    "Marker Swarm Plots" => "plot_tools/marker_plots.md",
    "Yield Strength Plots" => "plot_tools/yield_strength_plots.md",
    "Stokes Convergence Plots" => "plot_tools/stokes_convergence_plots.md"
]
data_types_pages = [
    "Data Types" => "data_types/data_types.md"
]
model_data_pages = [
    "ModelData" => "model_data/model_data.md",
    "Collections" => [
        "Benchmarks" => "model_data/collections/benchmarks.md",
        "BoundaryConditions" => "model_data/collections/bcs.md",
        "Carbonate" => "model_data/collections/carbonate.md",
        "Conversion" => "model_data/collections/conversion.md",
        "Geometry" => "model_data/collections/geometry.md",
        "Gravity" => "model_data/collections/gravity.md",
        "Grids" => "model_data/collections/grids.md",
        "HeatEquation" => "model_data/collections/heat_equation.md",
        "Interpolation" => "model_data/collections/interpolation.md",
        "Markers" => "model_data/collections/markers.md",
        "Materials" => "model_data/collections/materials.md",
        "Melting" => "model_data/collections/melting.md",
        "StokesContinuity" => "model_data/collections/stokes_continuity.md",
        "Timestep" => "model_data/collections/timestep.md",
        "Topography" => "model_data/collections/topography.md"
    ]
]

benchmark_pages = [
    "Benchmarks" => "benchmarks/benchmarks.md",
    "Single Execution" => "benchmarks/single_execution.md",
    "Batch Execution" => "benchmarks/batch_execution.md",
    "Results" => [
        "Couette Flow with Viscous Heating" => "benchmarks/models/couette_flow_viscous_heating.md",
        "Channel Flow Variable Conductivity" => "benchmarks/models/channel_flow_variable_conductivity.md",
        "Channel Flow Non Steady Temperature" => "benchmarks/models/channel_flow_non_steady_temperature.md",
        "Channel Flow Non Newtonian" => "benchmarks/models/channel_flow_non_newtonian.md",
        "Solid Body Rotation" => "benchmarks/models/solid_body_rotation.md",
        "Rayleigh Taylor Instability" => "benchmarks/models/rayleigh_taylor_instability.md",
        "Elastic Slab" => "benchmarks/models/elastic_slab.md",
        "Viscoelastic Stress Buildup" => "benchmarks/models/viscoelastic_stress_buildup.md",
        "Box Convection Isoviscous 1a" => "benchmarks/models/box_convection_isoviscous_1a.md",
        "Plasticity Benchmark Kaus10" => "benchmarks/models/plasticity_benchmark_kaus10.md",
        "Seafloor Spreading" => "benchmarks/models/seafloor_spreading.md",
        "Flexure Triangular Hole" => "benchmarks/models/flexure_triangular_hole.md",
        "Viscoelastic Extension" => "benchmarks/models/viscoelastic_extension.md",
        "Viscoelastic Contraction" => "benchmarks/models/viscoelastic_contraction.md",
        "Simple Sedimentation" => "benchmarks/models/simple_sedimentation.md"
    ]
]

units_pages = [
    "Units" => "units/units.md"
]

bibliography = "References" => "bibliography/bibliography.md"

PAGES = [
    introduction,
    "Installation" => install_pages,
    "Models" => model_pages,
    "Methods" => methods_pages,
    "Thermo-mechanical Algorithm" => thermo_mechanical_algorithm_pages,
    "EarthBox API" => earthbox_api_pages,
    "Material Library" => material_library_pages,
    "Plot Tools" => plot_tools_pages,
    "Data Types" => data_types_pages,
    "Model Data" => model_data_pages,
    "Units" => units_pages,
    "Benchmarks" => benchmark_pages,
    bibliography
]

makedocs(
    sitename = "EarthBox.jl", 
    remotes  = nothing,
    plugins  = [bib],
    pages    = PAGES,
    repo     = "https://github.com/eakneller/EarthBox.jl",
    #warnonly = [:cross_references],
    #draft = false,        # set to true to speed things up
    #doctest = true,       # set to false to speed things up
    #clean = true,
    format   = Documenter.HTML(
        collapselevel = 1, 
        assets = ["assets/custom.css", "assets/equation-refs.js"],
        search_size_threshold_warn = 1_000_000, # in bytes
    )
)

deploydocs(
    repo      = "github.com/eakneller/EarthBox.jl.git",
    devbranch = "main",               # your development branch
    # push_preview = true,            # optional: PR previews
    #forcepush = true,
)