module RightParts

import EarthBox.Arrays.ArrayTypes.ScalarArray3D: grid_array3D
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: grid_array2D
import ....LevelManager: LevelData, LevelData2d
import ....Domain
import ..ModelParametersManager: ModelParameters
import ..Density: average_density_basic_xz, average_density_basic_x

function build_right_parts_model(
    level_data::LevelData,
    rho::Array{Float64, 3},
    gravitational_acceleration::Float64
)::Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}}
    ynum = level_data.grid.parameters.geometry.ynum.value
    xnum = level_data.grid.parameters.geometry.xnum.value
    znum = level_data.grid.parameters.geometry.znum.value
    RX = grid_array3D(ynum, xnum, znum, Val(:vx))
    RY = grid_array3D(ynum, xnum, znum, Val(:vy))
    RZ = grid_array3D(ynum, xnum, znum, Val(:vz))
    for k = 1:znum
        for j = 1:xnum
            for i = 1:ynum
                if Domain.is_right_part_of_x_stokes_equation(i, j, k, ynum)
                    RX[i,j,k] = 0.0
                end
                if Domain.is_right_part_of_y_stokes_equation(i, j, k, ynum)
                    RY[i,j,k] = calculate_right_part_of_y_stokes_equation_3d(
                        i, j, k, rho, gravitational_acceleration)
                end
                if Domain.is_right_part_of_z_stokes_equation(i, j, k, ynum)
                    RZ[i,j,k] = 0.0
                end
            end
        end
    end
    return RX, RY, RZ
end

function build_right_parts_model(
    level_data::LevelData2d,
    rho::Array{Float64, 2},
    gravitational_acceleration::Float64
)::Tuple{Array{Float64, 2}, Array{Float64, 2}}
    ynum = level_data.grid.parameters.geometry.ynum.value
    xnum = level_data.grid.parameters.geometry.xnum.value
    RX = grid_array2D(ynum, xnum, Val(:vx))
    RY = grid_array2D(ynum, xnum, Val(:vy))
    for j = 1:xnum
        for i = 1:ynum
            if Domain.is_right_part_of_x_stokes_equation_2d(i, j, xnum)
                RX[i,j] = 0.0
            end
            if Domain.is_right_part_of_y_stokes_equation_2d(i, j, ynum)
                RY[i,j] = calculate_right_part_of_y_stokes_equation_2d(
                    i, j, rho, gravitational_acceleration)
            end
        end
    end
    return RX, RY
end

function calculate_right_part_of_y_stokes_equation_3d(
    i::Int64,
    j::Int64,
    k::Int64,
    rho::Array{Float64, 3},
    g::Float64
)::Float64
    return -g * average_density_basic_xz(i, j, k, rho)
end

function calculate_right_part_of_y_stokes_equation_2d(
    i::Int64,
    j::Int64,
    rho::Array{Float64, 2},
    g::Float64
)::Float64
    return -g * average_density_basic_x(i, j, rho)
end

end # module