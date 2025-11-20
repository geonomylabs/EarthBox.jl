module StokesSinker

include("model_data/ModelParametersManager.jl")
include("initial_conditions/Density.jl")
include("initial_conditions/Pressure.jl")
include("initial_conditions/Viscosity.jl")
include("initial_conditions/RightParts.jl")

import EarthBox.StaggeredGrid.Options: option_names
import EarthBox.Arrays.ArrayTypes.ScalarArray3D: grid_array3D
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: grid_array2D
import ...MultigridDataManager: MultigridData3d, MultigridData2d
import ...MultigridDataManager: MultigridDataManager
import ...ScalingFactors: calculate_residual_scaling_factors
import ...MultigridSolver: run_multigrid_solver
import .ModelParametersManager

function run_stokes_sinker(;
    make_plots::Bool=false, 
    use_multimulti::Bool=false, 
    model_type::Symbol=:ThreeDimensional
)::Nothing
    # Create the multigrid data structure
    if model_type == :ThreeDimensional
        grid = (
            ynum=97, xnum=97, znum=97,
            ysize=100000.0, xsize=100000.0, zsize=100000.0, 
            grid_type=option_names.UniformGrid
            )
        vcycle = (
            use_multimulti=use_multimulti,
            nvcycles=250, nvcycles_viscosity_jump=15, nviscosity_jumps=11.475031046555408,
            smoothing_iterations_on_finest_level=5, convergence_criterion=1e-20,
            make_plots=make_plots
            )
        relaxation = (
            relax_stokes=0.9, relax_velocity=1.0, relax_continuity=0.3,
            relax_pressure=1.0
            )
        multigrid_data = MultigridData3d(
            grid=grid, vcycle=vcycle, relaxation=relaxation, pressure_bc=0.0,
            gravitational_acceleration=9.81
            )
    elseif model_type == :TwoDimensional
        grid = (
            ynum=97, xnum=97,
            ysize=100000.0, xsize=100000.0,
            grid_type=option_names.UniformGrid
            )
        vcycle = (
            use_multimulti=use_multimulti,
            nvcycles=500, nvcycles_viscosity_jump=17, nviscosity_jumps=12.0, #11.475031046555408,
            smoothing_iterations_on_finest_level=5, convergence_criterion=1e-5,
            make_plots=make_plots
            )
        relaxation = (
            relax_stokes=0.9, relax_velocity=1.0, relax_continuity=0.3,
            relax_pressure=1.0
            )
        multigrid_data = MultigridData2d(
            grid=grid, vcycle=vcycle, relaxation=relaxation, pressure_bc=0.0,
            gravitational_acceleration=9.81
            )
    end
    # Set viscosity scaling parameters for multimulti. The default values above 
    # are for normal multigrid.
    if use_multimulti
        multigrid_data.viscosity_scaling.nvcycles_viscosity_jump = 3
        multigrid_data.viscosity_scaling.nviscosity_jumps = 5.9999
        multigrid_data.vcycle.convergence_criterion = 1e-10
    end

    # Define initial conditions, residual scaling factors and initial guesses
    model_parameters = ModelParametersManager.create_model_parameters(model_type=model_type)
    model_parameters.viscosity_block = 1e+26
    (
        multigrid_data.residual_scaling_factors.stokesscale,
        multigrid_data.residual_scaling_factors.continscale
    ) = calculate_residual_scaling_factors(
        model_parameters.density_block, model_parameters.viscosity_medium, 
        multigrid_data.gravitational_acceleration, 
        multigrid_data.level_vector[1].grid.parameters.geometry.xstpavg.value
        )
    if model_type == :ThreeDimensional
        # Initial conditions
        rho_ini = Density.build_density_model(
            multigrid_data.level_vector[1], model_parameters
            )
        etan_ini, etaxy_ini, etaxz_ini, etayz_ini = Viscosity.build_viscosity_model(
            multigrid_data.level_vector[1], model_parameters
            )
        RX_ini, RY_ini, RZ_ini = RightParts.build_right_parts_model(
            multigrid_data.level_vector[1], rho_ini, 
            multigrid_data.gravitational_acceleration
            )
        # Initial guesses
        pr_guess = Pressure.build_pressure_model(
            multigrid_data.level_vector[1], rho_ini, 
            multigrid_data.pressure_bc, multigrid_data.gravitational_acceleration
            )
        vx_guess, vy_guess, vz_guess = get_initial_velocity_guess(multigrid_data)

        MultigridDataManager.initialize_multigrid_data(
            multigrid_data, rho_ini, etan_ini, etaxy_ini, etaxz_ini, etayz_ini,
            vx_guess, vy_guess, vz_guess, pr_guess, RX_ini, RY_ini, RZ_ini
            )
    elseif model_type == :TwoDimensional
        # Initial conditions
        rho_ini = Density.build_density_model(
            multigrid_data.level_vector[1], model_parameters
            )
        etan_ini, etas_ini = Viscosity.build_viscosity_model(
            multigrid_data.level_vector[1], model_parameters
            )
        RX_ini, RY_ini = RightParts.build_right_parts_model(
            multigrid_data.level_vector[1], rho_ini, 
            multigrid_data.gravitational_acceleration
            )
        # Initial guesses
        pr_guess = Pressure.build_pressure_model(
            multigrid_data.level_vector[1], rho_ini, 
            multigrid_data.pressure_bc, multigrid_data.gravitational_acceleration
            )
        vx_guess, vy_guess = get_initial_velocity_guess(multigrid_data)

        MultigridDataManager.initialize_multigrid_data(
            multigrid_data, rho_ini, etan_ini, etas_ini,
            vx_guess, vy_guess, pr_guess, RX_ini, RY_ini
            )
    end

    run_multigrid_solver(multigrid_data)

    return nothing
end

function get_initial_velocity_guess(
    multigrid_data::MultigridData3d
)::Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}}
    xnum = multigrid_data.level_vector[1].grid.parameters.geometry.xnum.value
    ynum = multigrid_data.level_vector[1].grid.parameters.geometry.ynum.value
    znum = multigrid_data.level_vector[1].grid.parameters.geometry.znum.value
    vx_guess = grid_array3D(ynum, xnum, znum, Val(:vx))
    vy_guess = grid_array3D(ynum, xnum, znum, Val(:vy))
    vz_guess = grid_array3D(ynum, xnum, znum, Val(:vz))
    return vx_guess, vy_guess, vz_guess
end

function get_initial_velocity_guess(
    multigrid_data::MultigridData2d
)::Tuple{Array{Float64, 2}, Array{Float64, 2}}
    xnum = multigrid_data.level_vector[1].grid.parameters.geometry.xnum.value
    ynum = multigrid_data.level_vector[1].grid.parameters.geometry.ynum.value
    vx_guess = grid_array2D(ynum, xnum, Val(:vx))
    vy_guess = grid_array2D(ynum, xnum, Val(:vy))
    return vx_guess, vy_guess
end

end