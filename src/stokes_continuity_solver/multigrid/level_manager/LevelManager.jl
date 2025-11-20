module LevelManager

import EarthBox.Arrays.ArrayTypes.ScalarArray2D: ScalarArray2DState
import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.Arrays.ArrayTypes.ScalarArray3D: ScalarArray3DState
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import EarthBox.StaggeredGrid: InitManager
import ..GridMappingManager: MappingGroup, MappingGroup2d
import ..GridMappingManager: calculate_mapping_group!

const MAX_LEVELS = 10

mutable struct LevelData
    # level ID equal to multigrid level number with larger numbers for coarser grids
    level_id::Int64
    # velocity arrays m/s
    vx::ScalarArray3DState
    vy::ScalarArray3DState
    vz::ScalarArray3DState
    # pressure array Pa
    pr::ScalarArray3DState
    # density array kg/m^3
    rho::ScalarArray3DState
    # right-hand-side-part arrays on 3D grid
    RX::ScalarArray3DState
    RY::ScalarArray3DState
    RZ::ScalarArray3DState
    RC::ScalarArray3DState
    # viscosity arrays Pa.s
    etan::ScalarArray3DState
    etaxy::ScalarArray3DState
    etaxz::ScalarArray3DState
    etayz::ScalarArray3DState
    # backup viscosity arrays Pa.s
    etano::ScalarArray3DState
    etaxyo::ScalarArray3DState
    etaxzo::ScalarArray3DState
    etayzo::ScalarArray3DState
    # 3D grid data
    grid::Grids3d
    # fine to coarse grid mapping
    fine_to_coarse_mapping::MappingGroup
    # coarse to fine grid mapping
    coarse_to_fine_mapping::MappingGroup
    # grid type symbol
    grid_type::Symbol
end

"""
    LevelData(level_id, ynum, xnum, znum, ysize, xsize, zsize)

Construct a new LevelData with the specified parameters.

# Arguments
- `level_id::Int64`: Level ID
- `ynum::Int64`: Number of grid points in y-direction
- `xnum::Int64`: Number of grid points in x-direction
- `znum::Int64`: Number of grid points in z-direction
- `ysize::Float64`: Domain size in y-direction
- `xsize::Float64`: Domain size in x-direction
- `zsize::Float64`: Domain size in z-direction
- `grid_type::Symbol`: Grid type

# Returns
- `LevelData`: A new instance with initialized arrays and grid data for a given multigrid level
"""
function LevelData(
    level_id::Int64,
    ynum::Int64,
    xnum::Int64,
    znum::Int64,
    ysize::Float64,
    xsize::Float64,
    zsize::Float64,
    grid_type::Symbol
)::LevelData
    @assert ynum > 2 "Number of grid points in y-direction must be greater than 2"
    @assert xnum > 2 "Number of grid points in x-direction must be greater than 2"
    @assert znum > 2 "Number of grid points in z-direction must be greater than 2"
    @assert ysize > 0.0 "Domain size in y-direction must be positive"
    @assert xsize > 0.0 "Domain size in x-direction must be positive"
    @assert zsize > 0.0 "Domain size in z-direction must be positive"
    @assert level_id >= 0 "Level ID must be positive"

    grid = initialize_grid!(ynum, xnum, znum, ysize, xsize, zsize, grid_type)
    vx = ScalarArray3DState(ynum, xnum, znum, "vx", "m/s", "vx", "Velocity x-component on vx grid")
    vy = ScalarArray3DState(ynum, xnum, znum, "vy", "m/s", "vy", "Velocity y-component on vy grid")
    vz = ScalarArray3DState(ynum, xnum, znum, "vz", "m/s", "vz", "Velocity z-component on vz grid")
    pr = ScalarArray3DState(ynum, xnum, znum, "pr", "Pa", "pressure", "Pressure on pressure grid")
    
    rho = ScalarArray3DState(ynum, xnum, znum, "rho", "kg/m^3", "basic", "Density on basic grid")
    RX = ScalarArray3DState(ynum, xnum, znum, "RX", "None", "vx", "Right-hand-side-part on vx grid")
    RY = ScalarArray3DState(ynum, xnum, znum, "RY", "None", "vy", "Right-hand-side-part on vy grid")
    RZ = ScalarArray3DState(ynum, xnum, znum, "RZ", "None", "vz", "Right-hand-side-part on vz grid")
    RC = ScalarArray3DState(ynum, xnum, znum, "RC", "None", "pressure", "Right-hand-side-part on pressure grid")
    etan = ScalarArray3DState(ynum, xnum, znum, "etan", "Pa.s", "pressure", "Viscosity on pressure grid")
    etaxy = ScalarArray3DState(ynum, xnum, znum, "etaxy", "Pa.s", "shearxy", "Viscosity on shearxy grid")
    etaxz = ScalarArray3DState(ynum, xnum, znum, "etaxz", "Pa.s", "shearxz", "Viscosity on shearxz grid")
    etayz = ScalarArray3DState(ynum, xnum, znum, "etayz", "Pa.s", "shearyz", "Viscosity on shearyz grid")
    etano = ScalarArray3DState(ynum, xnum, znum, "etano", "Pa.s", "pressure", "Backup viscosity on pressure grid")
    etaxyo = ScalarArray3DState(ynum, xnum, znum, "etaxyo", "Pa.s", "shearxy", "Backup viscosity on shearxy grid")
    etaxzo = ScalarArray3DState(ynum, xnum, znum, "etaxzo", "Pa.s", "shearxz", "Backup viscosity on shearxz grid")
    etayzo = ScalarArray3DState(ynum, xnum, znum, "etayzo", "Pa.s", "shearyz", "Backup viscosity on shearyz grid")
    fine_to_coarse_mapping = MappingGroup(vx, vy, vz, pr, etaxy, etaxz, etayz)
    coarse_to_fine_mapping = MappingGroup(vx, vy, vz, pr, etaxy, etaxz, etayz)

    return LevelData(
        level_id,
        vx, vy, vz, pr, rho, 
        RX, RY, RZ, RC,
        etan, etaxy, etaxz, etayz, 
        etano, etaxyo, etaxzo, etayzo,
        grid, fine_to_coarse_mapping, 
        coarse_to_fine_mapping,
        grid_type
    )
end

mutable struct LevelData2d
    # level ID equal to multigrid level number with larger numbers for coarser grids
    level_id::Int64
    # velocity arrays m/s
    vx::ScalarArray2DState
    vy::ScalarArray2DState
    # pressure array Pa
    pr::ScalarArray2DState
    # density array kg/m^3
    rho::ScalarArray2DState
    # right-hand-side-part arrays on 3D grid
    RX::ScalarArray2DState
    RY::ScalarArray2DState
    RC::ScalarArray2DState
    # viscosity arrays Pa.s
    etan::ScalarArray2DState
    etas::ScalarArray2DState
    # backup viscosity arrays Pa.s
    etano::ScalarArray2DState
    etaso::ScalarArray2DState
    # 2D grid data
    grid::Grids
    # fine to coarse grid mapping
    fine_to_coarse_mapping::MappingGroup2d
    # coarse to fine grid mapping
    coarse_to_fine_mapping::MappingGroup2d
    # grid type symbol
    grid_type::Symbol
end

"""
    LevelData(level_id, ynum, xnum, znum, ysize, xsize, zsize)

Construct a new LevelData with the specified parameters.

# Arguments
- `level_id::Int64`: Level ID
- `ynum::Int64`: Number of grid points in y-direction
- `xnum::Int64`: Number of grid points in x-direction
- `znum::Int64`: Number of grid points in z-direction
- `ysize::Float64`: Domain size in y-direction
- `xsize::Float64`: Domain size in x-direction
- `zsize::Float64`: Domain size in z-direction

# Returns
- `LevelData`: A new instance with initialized arrays and grid data for a given multigrid level
"""
function LevelData2d(
    level_id::Int64,
    ynum::Int64,
    xnum::Int64,
    ysize::Float64,
    xsize::Float64,
    grid_type::Symbol
)::LevelData2d
    @assert ynum > 2 "Number of grid points in y-direction must be greater than 2"
    @assert xnum > 2 "Number of grid points in x-direction must be greater than 2"
    @assert ysize > 0.0 "Domain size in y-direction must be positive"
    @assert xsize > 0.0 "Domain size in x-direction must be positive"
    @assert level_id >= 0 "Level ID must be positive"

    grid = initialize_grid!(ynum, xnum, ysize, xsize, grid_type)
    vx = ScalarArray2DState(ynum, xnum, "vx", "m/s", "vx", "Velocity x-component on vx grid")
    vy = ScalarArray2DState(ynum, xnum, "vy", "m/s", "vy", "Velocity y-component on vy grid")
    pr = ScalarArray2DState(ynum, xnum, "pr", "Pa", "pressure", "Pressure on pressure grid")
    rho = ScalarArray2DState(ynum, xnum, "rho", "kg/m^3", "basic", "Density on basic grid")
    RX = ScalarArray2DState(ynum, xnum, "RX", "None", "vx", "Right-hand-side-part on vx grid")
    RY = ScalarArray2DState(ynum, xnum, "RY", "None", "vy", "Right-hand-side-part on vy grid")
    RC = ScalarArray2DState(ynum, xnum, "RC", "None", "pressure", "Right-hand-side-part on pressure grid")
    etan = ScalarArray2DState(ynum, xnum, "etan", "Pa.s", "pressure", "Viscosity on pressure grid")
    etas = ScalarArray2DState(ynum, xnum, "etas", "Pa.s", "basic", "Viscosity on basic grid")
    etano = ScalarArray2DState(ynum, xnum, "etano", "Pa.s", "pressure", "Backup viscosity on pressure grid")
    etaso = ScalarArray2DState(ynum, xnum, "etaso", "Pa.s", "basic", "Backup viscosity on basic grid")
    fine_to_coarse_mapping = MappingGroup2d(vx, vy, pr, etas)
    coarse_to_fine_mapping = MappingGroup2d(vx, vy, pr, etas)

    return LevelData2d(
        level_id,
        vx, vy, pr, rho, 
        RX, RY, RC, 
        etan, etas, 
        etano, etaso,
        grid, fine_to_coarse_mapping, 
        coarse_to_fine_mapping,
        grid_type
    )
end

function initialize_grid!(
    ynum::Int64,
    xnum::Int64,
    znum::Int64,
    ysize::Float64,
    xsize::Float64,
    zsize::Float64,
    option_name::Symbol
)::Grids3d
    grid = Grids3d(ynum=ynum, xnum=xnum, znum=znum, ysize=ysize, xsize=xsize, zsize=zsize)
    InitManager.initialize!(grid, Val(option_name))
    return grid
end

function initialize_grid!(
    ynum::Int64,
    xnum::Int64,
    ysize::Float64,
    xsize::Float64,
    option_name::Symbol
)::Grids
    grid = Grids(ynum, xnum, ysize, xsize)
    InitManager.initialize!(grid, Val(option_name))
    return grid
end

function save_original_viscosity!(level_data::LevelData)::Nothing
    level_data.etaxyo.array = deepcopy(level_data.etaxy.array)
    level_data.etaxzo.array = deepcopy(level_data.etaxz.array)
    level_data.etayzo.array = deepcopy(level_data.etayz.array)
    level_data.etano.array = deepcopy(level_data.etan.array)
    return nothing
end

function save_original_viscosity!(level_data::LevelData2d)::Nothing
    level_data.etaso.array = deepcopy(level_data.etas.array)
    level_data.etano.array = deepcopy(level_data.etan.array)
    return nothing
end

function save_original_viscosity!(level_vector::Vector{LevelData}, level_id::Int64)::Nothing
    level_vector[level_id].etaxyo.array = deepcopy(level_vector[level_id].etaxy.array)
    level_vector[level_id].etaxzo.array = deepcopy(level_vector[level_id].etaxz.array)
    level_vector[level_id].etayzo.array = deepcopy(level_vector[level_id].etayz.array)
    level_vector[level_id].etano.array = deepcopy(level_vector[level_id].etan.array)
    return nothing
end

function save_original_viscosity!(level_vector::Vector{LevelData2d}, level_id::Int64)::Nothing
    level_vector[level_id].etaso.array = deepcopy(level_vector[level_id].etas.array)
    level_vector[level_id].etano.array = deepcopy(level_vector[level_id].etan.array)
    return nothing
end


function calculate_multigrid_levels_3d(;
    ynum::Int64, 
    xnum::Int64, 
    znum::Int64,
)::Tuple{Int64, Vector{Int64}, Vector{Int64}, Vector{Int64}}
    # Initialize arrays to store dimensions at each level
    ynums = [ynum]
    xnums = [xnum] 
    znums = [znum]
    # Keep dividing dimensions by 2 until one reaches 4
    while minimum([ynums[end], xnums[end], znums[end]]) > 4
        ynum_new = ceil(Int, ynums[end]/2)
        xnum_new = ceil(Int, xnums[end]/2)
        znum_new = ceil(Int, znums[end]/2)
        min_num = min(ynum_new, xnum_new, znum_new)
        if min_num < 4 || length(ynums) >= MAX_LEVELS
            break
        end
        push!(ynums, ynum_new)
        push!(xnums, xnum_new)
        push!(znums, znum_new)
    end
    levelnum = length(ynums)
    return levelnum, ynums, xnums, znums
end

function calculate_multigrid_levels_2d(;
    ynum::Int64, 
    xnum::Int64, 
)::Tuple{Int64, Vector{Int64}, Vector{Int64}}
    # Initialize arrays to store dimensions at each level
    ynums = [ynum]
    xnums = [xnum] 
    # Keep dividing dimensions by 2 until one reaches 4
    while minimum([ynums[end], xnums[end]]) > 4
        ynum_new = ceil(Int, ynums[end]/2)
        xnum_new = ceil(Int, xnums[end]/2)
        min_num = min(ynum_new, xnum_new)
        if min_num < 4 || length(ynums) >= MAX_LEVELS
            break
        end
        push!(ynums, ynum_new)
        push!(xnums, xnum_new)
    end
    levelnum = length(ynums)
    return levelnum, ynums, xnums
end

function make_level_vector(;
    ynum::Int64,
    xnum::Int64,
    znum::Union{Int64, Nothing}=nothing,
    ysize::Float64,
    xsize::Float64,
    zsize::Union{Float64, Nothing}=nothing,
    option_name::Symbol,
)::Union{Vector{LevelData}, Vector{LevelData2d}}
    if znum === nothing || zsize === nothing
        return make_level_vector_2d(
            ynum=ynum, xnum=xnum, 
            ysize=ysize, xsize=xsize, 
            option_name=option_name
            )
    else
        return make_level_vector_3d(
            ynum=ynum, xnum=xnum, znum=znum, 
            ysize=ysize, xsize=xsize, zsize=zsize, 
            option_name=option_name
            )
    end
end

function make_level_vector_3d(;
    ynum::Int64,
    xnum::Int64,
    znum::Int64,
    ysize::Float64,
    xsize::Float64,
    zsize::Float64,
    option_name::Symbol
)::Vector{LevelData}
    (
        levelnum, ynums, xnums, znums
    ) = calculate_multigrid_levels_3d(ynum=ynum, xnum=xnum, znum=znum)
    level_vector = [
        LevelData(i, ynums[i], xnums[i], znums[i], ysize, xsize, zsize, option_name) 
        for i in 1:levelnum]
    make_fine_to_coarse_mapping!(level_vector)
    make_coarse_to_fine_mapping!(level_vector)
    return level_vector
end

function make_level_vector_2d(;
    ynum::Int64,
    xnum::Int64,
    ysize::Float64,
    xsize::Float64,
    option_name::Symbol
)::Vector{LevelData2d}
    (
        levelnum, ynums, xnums
    ) = calculate_multigrid_levels_2d(ynum=ynum, xnum=xnum)
    level_vector = [
        LevelData2d(i, ynums[i], xnums[i], ysize, xsize, option_name) 
        for i in 1:levelnum]
    make_fine_to_coarse_mapping!(level_vector)
    make_coarse_to_fine_mapping!(level_vector)
    return level_vector
end

function make_fine_to_coarse_mapping!(
    level_vector::Union{Vector{LevelData}, Vector{LevelData2d}}
)::Nothing
    levelnum = length(level_vector)
    for i in 1:levelnum-1
        gridf = level_vector[i].grid
        gridc = level_vector[i+1].grid
        fine_to_coarse_mapping = level_vector[i].fine_to_coarse_mapping
        calculate_mapping_group!(fine_to_coarse_mapping, gridf, gridc)
    end
    return nothing
end

function make_coarse_to_fine_mapping!(
    level_vector::Union{Vector{LevelData}, Vector{LevelData2d}}
)::Nothing
    levelnum = length(level_vector)
    for i in levelnum:-1:2
        gridc = level_vector[i].grid
        gridf = level_vector[i-1].grid
        coarse_to_fine_mapping = level_vector[i].coarse_to_fine_mapping
        calculate_mapping_group!(coarse_to_fine_mapping, gridc, gridf)
    end
    return nothing
end

function reset_solution_grids_to_zero!(level_data::LevelData)::Nothing
    level_data.vx.array .= 0.0
    level_data.vy.array .= 0.0
    level_data.vz.array .= 0.0
    level_data.pr.array .= 0.0
    return nothing
end

function reset_solution_grids_to_zero!(level_data::LevelData2d)::Nothing
    level_data.vx.array .= 0.0
    level_data.vy.array .= 0.0
    level_data.pr.array .= 0.0
    return nothing
end

end # module