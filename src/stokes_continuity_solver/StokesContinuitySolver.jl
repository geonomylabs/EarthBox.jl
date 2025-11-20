module StokesContinuitySolver

include("utils/StokesInterpolate.jl")
include("utils/GravityFuncs.jl")
include("utils/Viscoelastic.jl")
include("utils/GridStress.jl")
include("utils/strain_rate/GridStrainRate.jl")
include("utils/strain_rate/MarkerStrainRate.jl")
include("utils/strain_rate/MarkerStrainRateRatio.jl")
include("utils/MarkerPressure.jl")
include("subgrid_stress/SubGridStress.jl")
include("solver/SolverManager.jl")
include("velocity_type/VelocityType.jl")
include("multigrid/MultigridManager.jl")

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.PrintFuncs: @timeit_memit
import EarthBox: TimeStep
import EarthBox: Kinematics
import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer: load_parameters!
import EarthBox.Rheology.PlasticFailure: GridViscoplasticityToMarkers
import EarthBox.ConfigurationManager.SolverConfig: SolverConfigState
import EarthBox.DebugPrint: print_stokes_system_of_equations_stats
import EarthBox.DebugPrint: print_stokes_solution_vector_stats
import EarthBox.DebugPrint: print_stokes_solution_stats
import SparseArrays: SparseMatrixCSC
import .SolverManager: StokesResiduals
import .SolverManager: StokesNorms
import .SolverManager: SystemSolver
import .SolverManager: StokesViscoelasticTerms
import .SolverManager: StokesContinuityRhs
import .SolverManager: StokesBuildManager
import .VelocityType
import .VelocityType: option_names as velocity_type_names
import .VelocityType: make_velocity_type_string
import .GridStrainRate
import .MarkerStrainRate
import .MarkerStrainRateRatio
import .MarkerPressure
import .GridStress
import .StokesInterpolate
import .GravityFuncs
import .SubGridStress
import .SubGridStress: StressCorrection

const PDATA = get_eb_parameters()

const DEBUG = false

struct ValidInputNames
    gravity_x::Symbol
    gravity_y::Symbol
    turn_off_gravity_y::Symbol
    nsteps_turn_off_gravity::Symbol
    iuse_interface_stabilization::Symbol
end

"""
    initialize!(
        model::ModelData;
        velocity_type::Union{Int, String, Symbol, Nothing}=nothing,
        kwargs...
    )::Nothing

Initialize the Stokes-continuity solver.

# Arguments
- `model::`[`ModelData`](@ref ModelData): The model data container containing 
    model parameters and arrays.

# Arguments
- `velocity_type::Union{Int, String, Symbol, Nothing}`: Controls the type of 
    velocity calculation. See the **Velocity Types** section below for information 
    on available velocity types. If `velocity_type` is nothing the option id from 
    the model data container will be used to define velocity type. The velocity type 
    is stored in the model data container as an integer ID (`itype_velocity`) and 
    a corresponding string name (`stype_velocity`). The velocity type parameters can be 
    accessed from the model data container as follows:
    - `itype_velocity = model.stokes_continuity.parameters.velocity_calc_option.itype_velocity.value`
    - `stype_velocity = model.stokes_continuity.parameters.velocity_calc_option.stype_velocity.value`

# Keyword Arguments
- `gravity_x::Float64`:
    - $(PDATA.gravity_x.description)
- `gravity_y::Float64`:
    - $(PDATA.gravity_y.description)
- `turn_off_gravity_y::Int64`:
    - $(PDATA.turn_off_gravity_y.description)
- `nsteps_turn_off_gravity::Int64`:
    - $(PDATA.nsteps_turn_off_gravity.description)
- `iuse_interface_stabilization::Int64`:
    - $(PDATA.iuse_interface_stabilization.description)

# Returns
- `Nothing`

---
# Velocity Types
---
$(make_velocity_type_string())

"""
function initialize!(
    model::ModelData;
    velocity_type::Union{Int, String, Symbol, Nothing}=nothing,
    kwargs...
)
    VelocityType.initialize!(model, velocity_type=velocity_type)
    load_parameters!(model, fieldnames(ValidInputNames); kwargs...)
end

function solve_viscoelastic_stokes_continuity_equations!(
    model::ModelData,
    solver_config::SolverConfigState
)::Nothing
    Li, Lj, Lv = build_system_of_equations(model)
    if VelocityType.is_velocity_from_stokes_solver(model)
        S, Ls = solve_system_of_equations(model, Li, Lj, Lv, solver_config)
        SystemSolver.process_stokes_solution!(model, S)
        calculate_residuals!(model, Ls)
    else
        S = zeros(1)
        Ls = zeros(1, 1)
    end
    return nothing
end

function build_system_of_equations(
    model::ModelData
)::Tuple{Vector{Int64}, Vector{Int64}, Vector{Float64}}
    Li, Lj, Lv = stokes_solver_build_steps(model)
    return Li, Lj, Lv
end

""" Stokes-Continuity solver build steps.

Algorithm:
1. Calculate visco-elasto-plastic stress and viscosity terms including
arrays etas0, etan0, sxy0 and sxx0 where:
    etas0 = etas1*(1-viscoelastic_factor),
    etan0 = etan1*(1-viscoelastic_factor),
    sxy0 = sxy1*viscoelastic_factor,
    sxx0 = sxx1*viscoelastic_factor,
    viscoelastic_factor =
        viscoplastic_viscosity
            /(viscoplastic_viscosity + timestep*shear_modulus),
and etas1, etan1, sxy1, sxx1 are viscoplastic viscosity and stress arrays
interpolated from markers to grid nodes.

The viscoelastic_factor comes from inserting the time-dependent finite
difference visco-elastic formulation for visco-elastic stress
into the Stokes equations (Gerya, 2010, pg 179). This leads to new
coefficients and right hand-side terms in the discretized form of the
equations. For example, the arrays sxx0 and sxy0 are the old stresses
that become part of the right-hand side terms (see step b below). Also,
etas0 and etan0 are used in in the coefficients of the large matrix for
the discretized Stokes equation-continuity equations.

2. Update right-hand-side part arrays RX1, RY1 and RC1 for x-Stokes,
y-Stokes and continuity equations, respectively where:
    RX1 = -gx*density - delta(sxx0)/dx - delta(sxy0)/dy,
    RY1 = -gy*density + delta(sxx0)/dy - delta(sxy0)/dx,
    RC = 2*sin(dilatation_angle)*plastic_strain_rate.
The right hand-side part arrays are used to define right-hand side
values for x-Stokes, y-Stokes and continuity equations during the
building of the system of equations for each respective unknown.
Note how the old stress with viscoelastic terms are added to RX1
and RY1.

3. Build the system of equations by looping over basic nodes (solution nodes)
using viscoplastic viscosity arrays etan0 and etas0 that include the 
viscoelastic_factor (see step a), right-hand side part arrays RX1, RY1 and
RC1, arrays that define staggered grid spacing (xstp, ystp, xstpc, ystpc), 
velocity boundary condition arrays (btopx, btopy, bbottomx, bbottomy, bleftx,
brightx), pressure boundary condition parameters (pressure_bc_mode, pressure_bc),
internal boundary condition arrays (bintern_zone, bintern_velocity) and 
pressure scaling factor (pscale) used to condition the system of equations.
The the large-matrix of the system of equation is defined with
arrays Li, Lj, Lv that store the i-index (rows), j-index (columns)
and coefficients for non-zero values, respectively. The build step
also produces array R, which is the right-hand side array of the
system of equations.

** Indices of arrays are not shown to simplify the description of each
step.

Returns:
- Li: Large matrix row indices for non-zero matrix elements
- Lj: Large matrix column indices for non-zero matrix elements
- Lv: Non-zero matrix values
"""
function stokes_solver_build_steps(
    model::ModelData
)::Tuple{Vector{Int64}, Vector{Int64}, Vector{Float64}}
    StokesViscoelasticTerms.calc_viscoelastic_terms!(model)
    StokesContinuityRhs.viscoelastic_rhs!(model)
    SystemSolver.initialize_solver!(model)
    @timeit_memit "Finished looping over nodes to calculate non-zero coefficients and rhs vector" begin
        system_vectors = StokesBuildManager.build_system_of_equations(model)
    end

    if VelocityType.is_velocity_from_solid_body_rotation(model)
        Kinematics.solid_body_rotation!(model)
    end

    return system_vectors.Li, system_vectors.Lj, system_vectors.Lv
end

function solve_system_of_equations(
    model::ModelData,
    Li::Vector{Int64},
    Lj::Vector{Int64},
    Lv::Vector{Float64},
    solver_config::SolverConfigState
)::Tuple{Vector{Float64}, SparseMatrixCSC{Float64,Int64}}
    S, Ls = SystemSolver.solve_system(model, Li, Lj, Lv, solver_config)
    return S, Ls
end

""" Stokes-continuity solution processing.

This function performs the following tasks:
1. Calculate residual arrays for x-Stokes, y-Stokes and continuity equations
2. Calculate grid strain rate arrays and spin array
3. Update main model time step
4. Execute post-Stokes-solver steps
5. Calculate Stokes solution norms
6. Back up solution arrays
"""
function stokes_solution_processing!(model::ModelData, inside_flags::Vector{Int8})
    if VelocityType.is_velocity_from_stokes_solver(model)
        @timeit_memit "Finished calculating residuals for Stokes-continuity" begin
            StokesResiduals.compute_stokes_residuals!(model)
        end
        @timeit_memit "Finished computing nonlinear residuals for Stokes-continuity" begin
            StokesResiduals.compute_stokes_nonlinear_residuals!(model)
        end
    end
    GridStrainRate.update_grid_strain_rate_and_spin!(model)
    main_time_loop = model.timestep.parameters.main_time_loop
    iupdate_timestep = main_time_loop.iupdate_timestep.value
    if iupdate_timestep == 1
        TimeStep.adjust_time_step_using_displacement_limit!(model)
    end
    execute_post_stokes_solver_steps!(model, inside_flags)
    StokesNorms.calculate_stokes_norms!(model)
    SystemSolver.backup_stokes_arrays!(model)
end

function calculate_residuals!(model::ModelData, Ls::SparseMatrixCSC{Float64,Int64})::Nothing
    @timeit_memit "Finished calculating nonlinear residuals for Stokes-continuity" begin
        StokesResiduals.stokes_calc_nonlinear_system_residual!(model, Ls)
    end
    return nothing
end

""" Manage grid-to-marker viscoplastic viscosity interpolation.

If node-based plasticity is used, interpolate viscoplastic viscosity from grid 
to marker array and record the plastic failure state.
"""
function interpolate_grid_viscoplastic_viscosity_to_markers!(
    model::ModelData, 
    inside_flags::Vector{Int8}
)::Nothing
        if node_based_plasticity_is_used(model)
            @timeit_memit "Finished interpolating grid viscoplastic viscosity to markers" begin
                GridViscoplasticityToMarkers.calc_marker_viscoplastic_viscosity!(model, inside_flags)
            end
        end
    return nothing
end

function node_based_plasticity_is_used(model::ModelData)::Bool
    return model.stokes_continuity.parameters.picard.itype_global.value == 1
end

function execute_post_stokes_solver_steps!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    update_grid_stress!(model)
    calculate_grid_stress_change!(model)
    update_marker_strain_rate!(model, inside_flags)
    calculate_marker_strain_rate_ratio!(model, inside_flags)
    calculate_marker_pressure!(model, inside_flags)
    return nothing
end

""" Forecast visco-elasto-plastic grid stress arrays sxx2 and sxy2

Forecast visco-elasto-plastic grid stress arrays `sxx2` and `sxy2`
using current visco-plastic viscosity arrays `etan1` and `etas1`, 
deviatoric grid strain rate arrays `exx` and `exy`, old deviatoric 
stress arrays `sxx1` and `sxy1` interpolated from markers prior to 
Picard loop and current model time step `timestep`.

"""
function update_grid_stress!(model::ModelData)::Nothing
    @timeit_memit "Finished forecasting viscoelastic grid stress" begin
        GridStress.forecast_viscoelastic_grid_stress!(model)
    end
    return nothing
end

""" Calculate grid stress change arrays dsxx and dsxy

Calculate grid stress change arrays `dsxy` and `dsxx` using 
forecasted deviatoric stress `sxx2` and `sxy2` and old deviatoric 
stress  `sxx1` and `sxy1` interpolated from markers prior to Picard 
iteration.

"""
function calculate_grid_stress_change!(model::ModelData)::Nothing
    @timeit_memit "Finished calculating grid stress change" begin
        GridStress.calculate_deviatoric_grid_stress_change!(model)
    end
    return nothing
end

""" Update marker strain rate arrays marker_exx and marker_exy

Interpolate grid strain rate arrays `exx` and `exy` based on the 
updated Stokes-continuity solutions to marker arrays `marker_exx` and 
`marker_exy`, respectively.

"""
function update_marker_strain_rate!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    @timeit_memit "Finished calculating marker strain rate" begin
        MarkerStrainRate.calculate_marker_strain_rate!(model, inside_flags)
    end
    return nothing
end

""" Calculate marker strain rate ratio array marker_sr_ratio

Calculate ratio of nodal-stress-based strain rate to velocity-based 
strain rate for each marker using interpolated marker strain rate arrays 
`marker_exx` and `marker_exy`, forecasted grid stress `sxx2` and 
`sxy2` and grid stress changes `dsxx` and `dsxy`. Strain rate 
ratios are stored in the marker array `marker_sr_ratio`.

"""
function calculate_marker_strain_rate_ratio!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    @timeit_memit "Finished calculating marker strain rate ratio" begin
        MarkerStrainRateRatio.calculate_marker_strain_rate_ratio!(model, inside_flags)
    end
    return nothing
end

""" Calculate marker pressure using updated staggered pressure solution

Calculate marker pressure using the updated staggered pressure
solution array `pr1`.

"""
function calculate_marker_pressure!(model::ModelData, inside_flags::Vector{Int8})::Nothing
    @timeit_memit "Finished calculating marker pressure" begin
        MarkerPressure.calculate_marker_pressure!(model, inside_flags)
    end
    return nothing
end

""" Update marker stress for grid stress changes and sub-grid stress diffusion

Update marker stress for grid stress changes and sub-grid stress
diffusion. Subgrid diffusion is only applied if the diffusion 
coefficient is greater than zero.

"""
function update_marker_stress_for_grid_and_subgrid_changes!(
    model::ModelData,
    inside_flags::Vector{Int8}
)::Nothing
    subgrid_diffusion_params = model.markers.parameters.subgrid_diffusion
    subgrid_diff_coef_stress = subgrid_diffusion_params.subgrid_diff_coef_stress.value
    @timeit_memit "Finished updating marker stress for grid and subgrid changes" begin
        if subgrid_diff_coef_stress > 0
            SubGridStress.Update.update_subgrid_stress!(model, inside_flags)
        end
    end
    @timeit_memit "Finished correcting marker stress for remaining stress" begin
        StressCorrection.apply_marker_stress_correction_for_remaining_grid_stress!(
            model, inside_flags)
    end
    return nothing
end

function interpolate_staggered_grid_velocity_to_basic_grid!(model::ModelData)::Nothing
    @timeit_memit "Finished interpolating staggered grid velocity to basic grid" begin
        StokesInterpolate.interpolate_stag_velocity_to_basic_grid!(model)
    end
    return nothing
end

function turn_off_gravity!(model::ModelData)
    if model.gravity.parameters.turn_off_gravity_y.value == 1
        @timeit_memit "Finished turning off gravity" begin
            GravityFuncs.turn_off_gravity!(model)
        end
    end
end

function run_multigrid3d_test(;make_plots::Bool=false)
    Multigrid3dTest.run_multigrid_test(make_plots=make_plots)
end

end # module 