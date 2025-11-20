module Density

import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.Arrays.ArrayTypes.ScalarArray3D: grid_array3D
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: grid_array2D
import ....LevelManager: LevelData, LevelData2d
import ....NormalizeCoordinates: normalize_coordinates_for_basic_node
import ..ModelParametersManager: ModelParameters, in_block

function build_density_model(
    level_data::LevelData,
    model_parameters::ModelParameters
)::Array{Float64, 3}
    ynum = level_data.grid.parameters.geometry.ynum.value
    xnum = level_data.grid.parameters.geometry.xnum.value
    znum = level_data.grid.parameters.geometry.znum.value
    rho = grid_array3D(ynum, xnum, znum, Val(:basic))
    for k = 1:znum
        for j = 1:xnum
            for i = 1:ynum
                rho[i,j,k] = calculate_density(i, j, k, level_data.grid, model_parameters)
            end
        end
    end
    return rho
end

function build_density_model(
    level_data::LevelData2d,
    model_parameters::ModelParameters
)::Array{Float64, 2}
    ynum = level_data.grid.parameters.geometry.ynum.value
    xnum = level_data.grid.parameters.geometry.xnum.value
    rho = grid_array2D(ynum, xnum, Val(:basic))
    for j = 1:xnum
        for i = 1:ynum
            rho[i,j] = calculate_density(i, j, level_data.grid, model_parameters)
        end
    end
    return rho
end

function calculate_density(
    i::Int64,
    j::Int64,
    k::Int64,
    grid::Grids3d,
    model_parameters::ModelParameters
)::Float64
    xx, yy, zz = normalize_coordinates_for_basic_node(i, j, k, grid)
    if in_block(model_parameters.block_geometry, xx, yy, zz)
        return model_parameters.density_block
    else
        return model_parameters.density_medium
    end
end

function calculate_density(
    i::Int64,
    j::Int64,
    grid::Grids,
    model_parameters::ModelParameters
)::Float64
    xx, yy = normalize_coordinates_for_basic_node(i, j, grid)
    if in_block(model_parameters.block_geometry, xx, yy)
        return model_parameters.density_block
    else
        return model_parameters.density_medium
    end
end

function average_density_basic_xz(
    i::Int64,
    j::Int64,
    k::Int64,
    rho::Array{Float64, 3}
)::Float64
    return (rho[i,j,k] + rho[i,j-1,k] + rho[i,j,k-1] + rho[i,j-1,k-1]) / 4.0
end

function average_density_basic_x(
    i::Int64,
    j::Int64,
    rho::Array{Float64, 2}
)::Float64
    return (rho[i,j] + rho[i,j-1]) / 2.0
end

end # module