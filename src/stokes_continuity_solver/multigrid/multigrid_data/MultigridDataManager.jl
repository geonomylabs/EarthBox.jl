module MultigridDataManager

include("structures/MultigridStructures.jl")

import .MultigridStructures: Counters, VcycleParameters, RelaxationParameters, 
    ViscosityScalingParameters, MeanResiduals
import .MultigridStructures: ResidualScalingFactors
import ..LevelManager
import ..LevelManager: LevelData, LevelData2d

mutable struct MultigridData2d
    gravitational_acceleration::Float64
    pressure_bc::Float64
    counters::Counters
    vcycle::VcycleParameters
    relaxation::RelaxationParameters
    viscosity_scaling::ViscosityScalingParameters
    mean_residuals::MeanResiduals
    residual_scaling_factors::ResidualScalingFactors
    # Number of Gauss-Seidel smoothing iterations for each resolution level (old name = smoothing_iterations).
    smoothing_iterations::Vector{Int64}
    level_vector::Vector{LevelData2d}
    level0::LevelData2d
end

""" Constructor for 2D multigrid data.

# Arguments
- `grid::NamedTuple`: Grid parameters.
    - `ynum::Int64`: Number of grid points in the y-direction.
    - `xnum::Int64`: Number of grid points in the x-direction.
    - `ysize::Float64`: Grid spacing in the y-direction.
    - `xsize::Float64`: Grid spacing in the x-direction.
    - `grid_type::Symbol`: Grid type. Options: `:UniformGrid`, `:TtypeRefinedGrid`.
- `vcycle::NamedTuple`: V-cycle parameters.
    - `use_multimulti::Bool`: Boolean flag controlling whether multi-multigrid is used.
    - `nvcycles::Int64`: Number of v-cycles.
    - `nvcycles_viscosity_jump::Int64`: Number of v-cycles between viscosity jumps.
    - `nviscosity_jumps::Int64`: Number of viscosity jumps used to gradually increase viscosity contrast.
    - `smoothing_iterations_on_finest_level::Int64`: Number of Gauss-Seidel smoothing iterations for the finest level.
    - `convergence_criterion::Float64`: Convergence criterion.
    - `make_plots::Bool`: Whether to make plots.
- `relaxation::NamedTuple`: Relaxation parameters.
    - `relax_stokes::Float64`: Gauss-Seidel relaxation parameter for Stokes equation.
    - `relax_velocity::Float64`: Gauss-Seidel relaxation parameter for velocity corrections.
    - `relax_continuity::Float64`: Gauss-Seidel relaxation parameter for continuity equation.
    - `relax_pressure::Float64`: Gauss-Seidel relaxation parameter for pressure corrections.
- `pressure_bc::Float64`: Pressure boundary condition.
- `gravitational_acceleration::Float64`: Gravitational acceleration.

# Returns
- `MultigridData2d`: 2D multigrid data.
"""
function MultigridData2d(;
    grid::NamedTuple,
    vcycle::NamedTuple,
    relaxation::NamedTuple,
    pressure_bc::Float64,
    gravitational_acceleration::Float64 = 9.81,
)::MultigridData2d

    check_input_tuple_names2d(grid, vcycle, relaxation)
    
    level_vector = LevelManager.make_level_vector(
        ynum=grid.ynum, xnum=grid.xnum, ysize=grid.ysize, 
        xsize=grid.xsize, option_name=grid.grid_type
        )

    smoothing_iterations = zeros(Int64, length(level_vector))
    set_smoothing_iterations(
        smoothing_iterations, vcycle.smoothing_iterations_on_finest_level)

    return MultigridData2d(
        gravitational_acceleration,
        pressure_bc,
        Counters(),
        VcycleParameters(vcycle),
        RelaxationParameters(relaxation),
        ViscosityScalingParameters(vcycle),
        MeanResiduals(vcycle.nvcycles),
        ResidualScalingFactors(),
        smoothing_iterations,
        level_vector,
        make_level0(level_vector[1])
    )
end

mutable struct MultigridData3d
    gravitational_acceleration::Float64
    pressure_bc::Float64
    counters::Counters
    vcycle::VcycleParameters
    relaxation::RelaxationParameters
    viscosity_scaling::ViscosityScalingParameters
    mean_residuals::MeanResiduals
    residual_scaling_factors::ResidualScalingFactors
    # Number of Gauss-Seidel smoothing iterations for each resolution level (old name = smoothing_iterations).
    smoothing_iterations::Vector{Int64}
    level_vector::Vector{LevelData}
    level0::LevelData
end

""" Constructor for 3D multigrid data.

# Arguments
- `grid::NamedTuple`: Grid parameters.
    - `ynum::Int64`: Number of grid points in the y-direction.
    - `xnum::Int64`: Number of grid points in the x-direction.
    - `znum::Int64`: Number of grid points in the z-direction.
    - `ysize::Float64`: Grid spacing in the y-direction.
    - `xsize::Float64`: Grid spacing in the x-direction.
    - `zsize::Float64`: Grid spacing in the z-direction.
    - `grid_type::Symbol`: Grid type. Options: `:UniformGrid`, `:TtypeRefinedGrid`.
- `vcycle::NamedTuple`: V-cycle parameters.
    - `use_multimulti::Bool`: Boolean flag controlling whether multi-multigrid is used.
    - `nvcycles::Int64`: Number of v-cycles.
    - `nvcycles_viscosity_jump::Int64`: Number of v-cycles between viscosity jumps.
    - `nviscosity_jumps::Int64`: Number of viscosity jumps used to gradually increase viscosity contrast.
    - `smoothing_iterations_on_finest_level::Int64`: Number of Gauss-Seidel smoothing iterations for the finest level.
    - `convergence_criterion::Float64`: Convergence criterion.
    - `make_plots::Bool`: Whether to make plots.
- `relaxation::NamedTuple`: Relaxation parameters.
    - `relax_stokes::Float64`: Gauss-Seidel relaxation parameter for Stokes equation.
    - `relax_velocity::Float64`: Gauss-Seidel relaxation parameter for velocity corrections.
    - `relax_continuity::Float64`: Gauss-Seidel relaxation parameter for continuity equation.
    - `relax_pressure::Float64`: Gauss-Seidel relaxation parameter for pressure corrections.
- `pressure_bc::Float64`: Pressure boundary condition.
- `gravitational_acceleration::Float64`: Gravitational acceleration.

# Returns
- `MultigridData3d`: 3D multigrid data.
"""
function MultigridData3d(;
    grid::NamedTuple,
    vcycle::NamedTuple,
    relaxation::NamedTuple,
    pressure_bc::Float64,
    gravitational_acceleration::Float64 = 9.81,
)::MultigridData3d

    check_input_tuple_names3d(grid, vcycle, relaxation)

    level_vector = LevelManager.make_level_vector(
        ynum=grid.ynum, xnum=grid.xnum, znum=grid.znum, 
        ysize=grid.ysize, xsize=grid.xsize, zsize=grid.zsize, 
        option_name=grid.grid_type
        )
    
    smoothing_iterations = zeros(Int64, length(level_vector))
    set_smoothing_iterations(
        smoothing_iterations, vcycle.smoothing_iterations_on_finest_level)

    return MultigridData3d(
        gravitational_acceleration,
        pressure_bc,
        Counters(),
        VcycleParameters(vcycle),
        RelaxationParameters(relaxation),
        ViscosityScalingParameters(vcycle),
        MeanResiduals(vcycle.nvcycles),
        ResidualScalingFactors(),
        smoothing_iterations,
        level_vector,
        make_level0(level_vector[1])
    )
end

function make_level0(
    level1::LevelData
)::LevelData
    # Get data from level 1
    xnum = level1.grid.parameters.geometry.xnum.value
    ynum = level1.grid.parameters.geometry.ynum.value
    znum = level1.grid.parameters.geometry.znum.value
    xsize = level1.grid.parameters.geometry.xsize.value
    ysize = level1.grid.parameters.geometry.ysize.value
    zsize = level1.grid.parameters.geometry.zsize.value
    grid_type = level1.grid_type
    level0 = LevelData(0, ynum, xnum, znum, ysize, xsize, zsize, grid_type)
    return level0
end

function make_level0(
    level1::LevelData2d
)::LevelData2d
    xnum = level1.grid.parameters.geometry.xnum.value
    ynum = level1.grid.parameters.geometry.ynum.value
    xsize = level1.grid.parameters.geometry.xsize.value
    ysize = level1.grid.parameters.geometry.ysize.value
    grid_type = level1.grid_type
    level0 = LevelData2d(0, ynum, xnum, ysize, xsize, grid_type)
    return level0
end

function copy_level1_density_and_viscosity_to_level0(
    multigrid_data::MultigridData3d
)::Nothing
    multigrid_data.level0.rho.array .= multigrid_data.level_vector[1].rho.array
    multigrid_data.level0.etan.array .= multigrid_data.level_vector[1].etan.array
    multigrid_data.level0.etaxy.array .= multigrid_data.level_vector[1].etaxy.array
    multigrid_data.level0.etaxz.array .= multigrid_data.level_vector[1].etaxz.array
    multigrid_data.level0.etayz.array .= multigrid_data.level_vector[1].etayz.array
    return level0
end

function copy_level1_density_and_viscosity_to_level0(
    multigrid_data::MultigridData2d
)::Nothing
    multigrid_data.level0.rho.array .= multigrid_data.level_vector[1].rho.array
    multigrid_data.level0.etan.array .= multigrid_data.level_vector[1].etan.array
    multigrid_data.level0.etas.array .= multigrid_data.level_vector[1].etas.array
    return nothing
end

function set_smoothing_iterations(
    smoothing_iterations::Vector{Int64},
    smoothing_iterations_on_finest_level::Int64,
)::Nothing
    smoothing_iterations[1] = smoothing_iterations_on_finest_level
    max_levels = length(smoothing_iterations)
    for i = 2:max_levels
        smoothing_iterations[i] = smoothing_iterations[i-1] * 2
    end
    return nothing
end

function initialize_multigrid_data(
    multigrid_data::MultigridData2d,
    rho_ini::Array{Float64, 2},
    etan_ini::Array{Float64, 2},
    etas_ini::Array{Float64, 2},
    vx_guess::Array{Float64, 2},
    vy_guess::Array{Float64, 2},
    pr_guess::Array{Float64, 2},
    RX_ini::Array{Float64, 2},
    RY_ini::Array{Float64, 2}
)::Nothing
    initialize_density_and_viscosity!(
        multigrid_data, rho_ini, etan_ini, etas_ini)
    initialize_pressure_and_rhs_parts!(multigrid_data, pr_guess, RX_ini, RY_ini)
    initialize_velocity!(multigrid_data, vx_guess, vy_guess)
    return nothing
end

function initialize_multigrid_data(
    multigrid_data::MultigridData3d,
    rho_ini::Array{Float64, 3},
    etan_ini::Array{Float64, 3},
    etaxy_ini::Array{Float64, 3},
    etaxz_ini::Array{Float64, 3},
    etayz_ini::Array{Float64, 3},
    vx_guess::Array{Float64, 3},
    vy_guess::Array{Float64, 3},
    vz_guess::Array{Float64, 3},
    pr_guess::Array{Float64, 3},
    RX_ini::Array{Float64, 3},
    RY_ini::Array{Float64, 3},
    RZ_ini::Array{Float64, 3}
)::Nothing
    initialize_density_and_viscosity!(
        multigrid_data, rho_ini, etan_ini, etaxy_ini, etaxz_ini, etayz_ini)
    initialize_pressure_and_rhs_parts!(multigrid_data, pr_guess, RX_ini, RY_ini, RZ_ini)
    initialize_velocity!(multigrid_data, vx_guess, vy_guess, vz_guess)
    return nothing
end

function initialize_density_and_viscosity!(
    multigrid_data::MultigridData3d,
    rho::Array{Float64, 3},
    etan::Array{Float64, 3},
    etaxy::Array{Float64, 3},
    etaxz::Array{Float64, 3},
    etayz::Array{Float64, 3}
)::Nothing
    multigrid_data.level_vector[1].rho.array .= rho
    multigrid_data.level0.rho.array .= rho
    multigrid_data.level_vector[1].etan.array .= etan
    multigrid_data.level_vector[1].etaxy.array .= etaxy
    multigrid_data.level_vector[1].etaxz.array .= etaxz
    multigrid_data.level_vector[1].etayz.array .= etayz
    multigrid_data.level0.etan.array .= etan
    multigrid_data.level0.etaxy.array .= etaxy
    multigrid_data.level0.etaxz.array .= etaxz
    multigrid_data.level0.etayz.array .= etayz
    set_initial_viscosity_scaling_parameters(multigrid_data)
    return nothing
end

function calculate_minimum_and_maximum_viscosity(
    multigrid_data::MultigridData3d
)::Tuple{Float64, Float64}
    etan_min = minimum(multigrid_data.level_vector[1].etan.array)
    etan_max = maximum(multigrid_data.level_vector[1].etan.array)
    etaxy_min = minimum(multigrid_data.level_vector[1].etaxy.array)
    etaxy_max = maximum(multigrid_data.level_vector[1].etaxy.array)
    etaxz_min = minimum(multigrid_data.level_vector[1].etaxz.array)
    etaxz_max = maximum(multigrid_data.level_vector[1].etaxz.array)
    etayz_min = minimum(multigrid_data.level_vector[1].etayz.array)
    etayz_max = maximum(multigrid_data.level_vector[1].etayz.array)
    eta_min = min(etan_min, etaxy_min, etaxz_min, etayz_min)
    eta_max = max(etan_max, etaxy_max, etaxz_max, etayz_max)
    return eta_min, eta_max
end

function calculate_minimum_and_maximum_viscosity(
    multigrid_data::MultigridData2d
)::Tuple{Float64, Float64}
    etan_min = minimum(multigrid_data.level_vector[1].etan.array)
    etan_max = maximum(multigrid_data.level_vector[1].etan.array)
    etas_min = minimum(multigrid_data.level_vector[1].etas.array)
    etas_max = maximum(multigrid_data.level_vector[1].etas.array)
    eta_min = min(etan_min, etas_min)
    eta_max = max(etan_max, etas_max)
    return eta_min, eta_max
end

function initialize_density_and_viscosity!(
    multigrid_data::MultigridData2d,
    rho::Array{Float64, 2},
    etan::Array{Float64, 2},
    etas::Array{Float64, 2}
)::Nothing
    multigrid_data.level_vector[1].rho.array .= rho
    multigrid_data.level0.rho.array .= rho
    multigrid_data.level_vector[1].etan.array .= etan
    multigrid_data.level_vector[1].etas.array .= etas
    multigrid_data.level0.etan.array .= etan
    multigrid_data.level0.etas.array .= etas
    set_initial_viscosity_scaling_parameters(multigrid_data)
    return nothing
end

function set_initial_viscosity_scaling_parameters(
    multigrid_data::Union{MultigridData3d, MultigridData2d}
)::Nothing
    viscosity_min, viscosity_max = calculate_minimum_and_maximum_viscosity(multigrid_data)
    println("Current viscosity range: $viscosity_min, $viscosity_max")
    multigrid_data.viscosity_scaling.viscosity_min_o = viscosity_min
    multigrid_data.viscosity_scaling.viscosity_max_o = viscosity_max
    multigrid_data.viscosity_scaling.viscosity_min_cur = viscosity_min
    multigrid_data.viscosity_scaling.viscosity_max_cur = viscosity_min
    multigrid_data.viscosity_scaling.viscosity_max_start = viscosity_min
    viscosity_max_factor = calculate_viscosity_scaling_factor(
        viscosity_min, viscosity_max, multigrid_data.viscosity_scaling.nviscosity_jumps)
    println("Calculated viscosity max factor: $viscosity_max_factor")
    println("N viscosity jumps: $(multigrid_data.viscosity_scaling.nviscosity_jumps)")
    multigrid_data.viscosity_scaling.viscosity_max_factor = viscosity_max_factor
    return nothing
end

function calculate_viscosity_scaling_factor(
    viscosity_min::Float64, 
    viscosity_max::Float64,
    nviscosity_jumps::Float64
)::Float64
    return 10^(log10(viscosity_max/viscosity_min)/nviscosity_jumps)
end

function initialize_pressure_and_rhs_parts!(
    multigrid_data::MultigridData3d,
    pr::Array{Float64, 3},
    RX::Array{Float64, 3},
    RY::Array{Float64, 3},
    RZ::Array{Float64, 3}
)::Nothing
    use_multimulti = multigrid_data.vcycle.use_multimulti
    if use_multimulti
        multigrid_data.level0.pr.array .= pr
        multigrid_data.level0.RX.array .= RX
        multigrid_data.level0.RY.array .= RY
        multigrid_data.level0.RZ.array .= RZ
    else
        multigrid_data.level_vector[1].pr.array .= pr
        multigrid_data.level_vector[1].RX.array .= RX
        multigrid_data.level_vector[1].RY.array .= RY
        multigrid_data.level_vector[1].RZ.array .= RZ
    end
    return nothing
end

function initialize_pressure_and_rhs_parts!(
    multigrid_data::MultigridData2d,
    pr::Array{Float64, 2},
    RX::Array{Float64, 2},
    RY::Array{Float64, 2}
)::Nothing
    use_multimulti = multigrid_data.vcycle.use_multimulti
    if use_multimulti
        multigrid_data.level0.pr.array .= pr
        multigrid_data.level0.RX.array .= RX
        multigrid_data.level0.RY.array .= RY
    else
        multigrid_data.level_vector[1].pr.array .= pr
        multigrid_data.level_vector[1].RX.array .= RX
        multigrid_data.level_vector[1].RY.array .= RY
    end
    return nothing
end

function initialize_velocity!(
    multigrid_data::MultigridData3d,
    vx::Array{Float64, 3},
    vy::Array{Float64, 3},
    vz::Array{Float64, 3}
)::Nothing
    multigrid_data.level_vector[1].vx.array .= vx
    multigrid_data.level_vector[1].vy.array .= vy
    multigrid_data.level_vector[1].vz.array .= vz
    # Do not initialize level zero since if level 0 is used level 1 acts as a correction
    # Validate this
    return nothing
end

function initialize_velocity!(
    multigrid_data::MultigridData2d,
    vx::Array{Float64, 2},
    vy::Array{Float64, 2}
)::Nothing
    multigrid_data.level_vector[1].vx.array .= vx
    multigrid_data.level_vector[1].vy.array .= vy
    # Do not initialize level zero since if level 0 is used level 1 acts as a correction
    # Validate this
    return nothing
end

function check_input_tuple_names2d(
    grid::NamedTuple,
    vcycle::NamedTuple,
    relaxation::NamedTuple,
)::Nothing
    check_grid_names2d(grid)
    check_vcycle_names(vcycle)
    check_relaxation_names(relaxation)
    return nothing
end

function check_input_tuple_names3d(
    grid::NamedTuple,
    vcycle::NamedTuple,
    relaxation::NamedTuple,
)::Nothing
    check_grid_names3d(grid)
    check_vcycle_names(vcycle)
    check_relaxation_names(relaxation)
    return nothing
end

function check_grid_names2d(grid::NamedTuple)::Nothing
    required_names = (:ynum, :xnum, :ysize, :xsize, :grid_type)
    check_required_names_for_named_tuple(grid, required_names)
    return nothing
end

function check_grid_names3d(grid::NamedTuple)::Nothing
    required_names = (:ynum, :xnum, :znum, :ysize, :xsize, :zsize, :grid_type)
    check_required_names_for_named_tuple(grid, required_names)
    return nothing
end

function check_vcycle_names(vcycle::NamedTuple)::Nothing
    required_names = (
        :use_multimulti, :nvcycles, :nvcycles_viscosity_jump, :nviscosity_jumps, 
        :smoothing_iterations_on_finest_level, :convergence_criterion, :make_plots
        )
    check_required_names_for_named_tuple(vcycle, required_names)
    return nothing
end

function check_relaxation_names(relaxation::NamedTuple)::Nothing
    required_names = (:relax_stokes, :relax_velocity, :relax_continuity, :relax_pressure)
    check_required_names_for_named_tuple(relaxation, required_names)
    return nothing
end

function check_required_names_for_named_tuple(
    input::NamedTuple,
    required_names::Tuple
)::Bool
    missing = filter(name -> !(name in keys(input)), required_names)
    if !isempty(missing)
        error("Input tuple is missing the required fields: $(missing). Check the input tuple names.")
    end
    return true
end

end # module