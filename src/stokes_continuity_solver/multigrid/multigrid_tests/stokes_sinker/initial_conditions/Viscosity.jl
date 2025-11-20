module Viscosity

import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.Arrays.ArrayTypes.ScalarArray3D: grid_array3D
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: grid_array2D
import ....LevelManager: LevelData, LevelData2d
import ....Domain
import ....NormalizeCoordinates
import ..ModelParametersManager: ModelParameters, in_block

function build_viscosity_model(
    level_data::LevelData,
    model_parameters::ModelParameters
)::Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}}
    ynum = level_data.grid.parameters.geometry.ynum.value
    xnum = level_data.grid.parameters.geometry.xnum.value
    znum = level_data.grid.parameters.geometry.znum.value
    etan = grid_array3D(ynum, xnum, znum, Val(:pressure))
    etaxy = grid_array3D(ynum, xnum, znum, Val(:shearxy))
    etaxz = grid_array3D(ynum, xnum, znum, Val(:shearxz))
    etayz = grid_array3D(ynum, xnum, znum, Val(:shearyz))
    for k = 1:znum
        for j = 1:xnum
            for i = 1:ynum
                if Domain.is_pressure_node(i, j, k, ynum, xnum, znum)
                    etan[i,j,k] = calculate_normal_viscosity(
                        i, j, k, level_data.grid, model_parameters)
                end
                if Domain.is_shearxy_node(k, znum)
                    etaxy[i,j,k] = calculate_shearxy_viscosity(
                        i, j, k, level_data.grid, model_parameters)
                end
                if Domain.is_shearxz_node(i, ynum)
                    etaxz[i,j,k] = calculate_shearxz_viscosity(
                        i, j, k, level_data.grid, model_parameters)
                end
                if Domain.is_shearyz_node(j, xnum)
                    etayz[i,j,k] = calculate_shearyz_viscosity(
                        i, j, k, level_data.grid, model_parameters)
                end
            end
        end
    end
    return etan, etaxy, etaxz, etayz
end

function build_viscosity_model(
    level_data::LevelData2d,
    model_parameters::ModelParameters
)::Tuple{Array{Float64, 2}, Array{Float64, 2}}
    ynum = level_data.grid.parameters.geometry.ynum.value
    xnum = level_data.grid.parameters.geometry.xnum.value
    etan = grid_array2D(ynum, xnum, Val(:pressure))
    etas = grid_array2D(ynum, xnum, Val(:basic))
    for j = 1:xnum
        for i = 1:ynum
            etas[i,j] = calculate_basic_shear_viscosity(
                i, j, level_data.grid, model_parameters)
            if Domain.is_pressure_node_2d(i, j, ynum, xnum)
                etan[i,j] = calculate_normal_viscosity(
                    i, j, level_data.grid, model_parameters)
            end
        end
    end
    return etan, etas
end

function calculate_normal_viscosity(
    i::Int64,
    j::Int64,
    k::Int64,
    grid::Grids3d,
    model_parameters::ModelParameters
)::Float64
    xx, yy, zz = NormalizeCoordinates.normalize_coordinates_for_pressure_node(i, j, k, grid)
    if in_block(model_parameters.block_geometry, xx, yy, zz)
        return model_parameters.viscosity_block
    else
        return model_parameters.viscosity_medium
    end
end

function calculate_normal_viscosity(
    i::Int64,
    j::Int64,
    grid::Grids,
    model_parameters::ModelParameters
)::Float64
    xx, yy = NormalizeCoordinates.normalize_coordinates_for_pressure_node(i, j, grid)
    if in_block(model_parameters.block_geometry, xx, yy)
        return model_parameters.viscosity_block
    else
        return model_parameters.viscosity_medium
    end
end

function calculate_basic_shear_viscosity(
    i::Int64,
    j::Int64,
    grid::Grids,
    model_parameters::ModelParameters
)::Float64
    xx, yy = NormalizeCoordinates.normalize_coordinates_for_basic_node(i, j, grid)
    if in_block(model_parameters.block_geometry, xx, yy)
        return model_parameters.viscosity_block
    else
        return model_parameters.viscosity_medium
    end
end

function calculate_shearxy_viscosity(
    i::Int64,
    j::Int64,
    k::Int64,
    grid::Grids3d,
    model_parameters::ModelParameters
)::Float64
    xx, yy, zz = NormalizeCoordinates.normalize_coordinates_for_shearxy_node(i, j, k, grid)
    if in_block(model_parameters.block_geometry, xx, yy, zz)
        return model_parameters.viscosity_block
    else
        return model_parameters.viscosity_medium
    end
end

function calculate_shearxz_viscosity(
    i::Int64,
    j::Int64,
    k::Int64,
    grid::Grids3d,
    model_parameters::ModelParameters
)::Float64
    xx, yy, zz = NormalizeCoordinates.normalize_coordinates_for_shearxz_node(i, j, k, grid)
    if in_block(model_parameters.block_geometry, xx, yy, zz)
        return model_parameters.viscosity_block
    else
        return model_parameters.viscosity_medium
    end
end

function calculate_shearyz_viscosity(
    i::Int64,
    j::Int64,
    k::Int64,
    grid::Grids3d,
    model_parameters::ModelParameters
)::Float64
    xx, yy, zz = NormalizeCoordinates.normalize_coordinates_for_shearyz_node(i, j, k, grid)
    if in_block(model_parameters.block_geometry, xx, yy, zz)
        return model_parameters.viscosity_block
    else
        return model_parameters.viscosity_medium
    end
end

end # module