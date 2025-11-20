module Pressure

import EarthBox.Arrays.ArrayTypes.ScalarArray3D: grid_array3D
import EarthBox.Arrays.ArrayTypes.ScalarArray2D: grid_array2D
import ....LevelManager: LevelData, LevelData2d
import ....Domain
import ..Density: average_density_basic_xz, average_density_basic_x
import ..ModelParametersManager: ModelParameters

function build_pressure_model(
    level_data::LevelData,
    rho::Array{Float64, 3},
    pressure_bc::Float64,
    gravitational_acceleration::Float64
)::Array{Float64, 3}
    ystp_pr = level_data.grid.arrays.pressure.ystp_pr.array
    ynum = level_data.grid.parameters.geometry.ynum.value
    xnum = level_data.grid.parameters.geometry.xnum.value
    znum = level_data.grid.parameters.geometry.znum.value
    pr = grid_array3D(ynum, xnum, znum, Val(:pressure))
    for k = 1:znum
        for j = 1:xnum
            for i = 1:ynum
                if Domain.is_associated_with_internal_pressure_node(i, j, k, ynum)
                    if i == 1
                        pr[i, j-1, k-1] = pressure_bc
                    else
                        # Calculate pressure change using average density in xz plane of basic grid
                        density_average = average_density_basic_xz(i, j, k, rho)
                        pr[i, j-1, k-1] = (
                              pr[i-1, j-1, k-1] 
                            + ystp_pr[i-1] 
                            * gravitational_acceleration 
                            * density_average
                            )
                    end
                end
            end
        end
    end
    return pr
end

function build_pressure_model(
    level_data::LevelData2d,
    rho::Array{Float64, 2},
    pressure_bc::Float64,
    gravitational_acceleration::Float64
)::Array{Float64, 2}
    ystp_pr = level_data.grid.arrays.pressure.ystp_pr.array
    ynum = level_data.grid.parameters.geometry.ynum.value
    xnum = level_data.grid.parameters.geometry.xnum.value
    pr = grid_array2D(ynum, xnum, Val(:pressure))
    for j = 1:xnum
        for i = 1:ynum
            if Domain.is_associated_with_internal_pressure_node_2d(i, j, ynum)
                if i == 1
                    pr[i, j-1] = pressure_bc
                else
                    # Calculate pressure change using average density in x direction basic grid
                    density_average = average_density_basic_x(i, j, rho)
                    pr[i, j-1] = (
                          pr[i-1, j-1] 
                        + ystp_pr[i-1] 
                        * gravitational_acceleration 
                        * density_average
                        )
                end
            end
        end
    end
    return pr
end

end