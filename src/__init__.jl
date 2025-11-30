# General utils - PrintFuncs must come first as it's used by MPIManager
include("utils/SysTools.jl")
include("utils/PrintFuncs.jl")
include("utils/PlotSettingsManager.jl")
include("utils/TupleTools.jl")
# MPI manager
include("mpi_manager/MPIManager.jl")
# Unit Conversion
include("utils/ConversionFuncs.jl")
# EarthBox Data Types
include("data_types/EarthBoxDtypes.jl")
include("data_types/colors/MaterialColorsContainer.jl")
include("data_types/arrays/Arrays.jl")
include("data_types/parameters/Parameters.jl")
include("data_types/parameters/SetParameters.jl")
include("data_types/parameters/ParameterGroupTools.jl")
include("data_types/parameters/ParameterRegistry.jl")
# Configuration
include("config/ConfigurationManager.jl")
# Material Library
include("material_library/MaterialLibraryCollection.jl")
# Input tools used read and validate input files
include("input_tools/InputTools.jl")
# General utils
include("utils/EBCopy.jl")
include("utils/DataStructures.jl")
include("utils/JLDTools.jl")
include("utils/FileMover.jl")
include("utils/Domain.jl")
include("utils/MathTools.jl")
include("utils/MarkerSwarm.jl")
include("utils/BuildSysTools.jl")
include("utils/GetArgs.jl")
include("utils/GetPaths.jl")
include("utils/GetModels.jl")
include("utils/UnitConversion.jl")
include("utils/EarthBoxPaths.jl")
include("utils/PathValidation.jl")
include("utils/SecurityUtils.jl")
include("utils/DictUtils.jl")
include("utils/FormatValueStrings.jl")
include("utils/UseOptionTools.jl")
include("utils/ViscoelasticFactor.jl")
include("utils/LoopInputStruct.jl")
include("utils/SedimentWaterInterface.jl")
include("utils/FindShallowest.jl")
include("utils/OptionTools.jl")
include("utils/RunTools.jl")
include("utils/TtypeCalculator.jl")
# Solidus and liquidus
include("solidus_liquidus/SolidusLiquidus.jl")
# EarthBox model data
include("data_types/parameters/ModelParametersInfo.jl")
include("model_data/ModelDataContainer.jl")
# EarthBox utils dependent on ModelDataContainer
include("utils/DebugPrint.jl")
include("utils/DebugPlots.jl")
include("utils/MobileWall.jl")
include("utils/InitializationTools.jl")
include("utils/GridFuncs.jl")
include("utils/MagmaFlushState.jl")
# Model case management tools
include("case_input_tools/CaseInputTools.jl")
# Structure
include("structure/ModelStructureManager.jl")
# Sediment thickness
include("utils/SedimentThickness.jl") # Dependent on ModelStructureManager and ModelDataContainer
# Compaction
include("compaction/Compaction.jl")
# Parallel solver
include("parallel_solver/ParallelSolver.jl")
# EarthBox model API
include("kinematics/Kinematics.jl")
include("staggered_grid/StaggeredGrid.jl")
include("interpolation/Interpolation.jl")
include("rheology/Rheology.jl")
include("markers/Markers.jl")
include("regrid/Regrid.jl")
include("advection/Advection.jl")
include("rock_properties/RockProperties.jl")
include("boundary_conditions/BoundaryConditions.jl")
include("time_step/TimeStep.jl") # Depends on BoundaryConditions
include("surface_processes/SurfaceProcesses.jl") # Depends on TimeStep: TimeStepCalc
include("stokes_continuity_solver/StokesContinuitySolver.jl")
include("global_plasticity_loop/GlobalPlasticityLoop.jl")
include("heat_solver/HeatSolver.jl")
include("melt_model/MeltModel.jl")
include("hydrothermal_circulation/HydrothermalCirculation.jl")
include("serpentinization/Serpentinization.jl")
include("material_geometry/MaterialGeometry.jl")
include("model_manager/ModelManager.jl")
include("time_loop/TimeLoop.jl")
include("potential_fields/gravity/Gravity.jl")
include("plot_tools/PlotToolsManager.jl")
include("benchmarks/BenchmarksManager.jl")
# Test tools
include("utils/TestTools.jl")
# Integrated tests
include("tests/TestManager.jl")