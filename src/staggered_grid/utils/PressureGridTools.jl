module PressureGridTools

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d

function calculate_grid_spacing!(grids::Grids)
    xnum = grids.parameters.geometry.xnum.value
    ynum = grids.parameters.geometry.ynum.value
    gridx_b = grids.arrays.basic.gridx_b.array
    gridy_b = grids.arrays.basic.gridy_b.array
    xstp_pr = grids.arrays.pressure.xstp_pr.array
    ystp_pr = grids.arrays.pressure.ystp_pr.array
    spacing_loop!(xnum, gridx_b, xstp_pr)
    spacing_loop!(ynum, gridy_b, ystp_pr)
    return nothing
end

function calculate_grid_spacing!(grids3d::Grids3d)::Nothing
    xnum = grids3d.parameters.geometry.xnum.value
    ynum = grids3d.parameters.geometry.ynum.value
    znum = grids3d.parameters.geometry.znum.value
    gridx_b = grids3d.arrays.basic.gridx_b.array
    gridy_b = grids3d.arrays.basic.gridy_b.array
    gridz_b = grids3d.arrays.basic.gridz_b.array
    xstp_pr = grids3d.arrays.pressure.xstp_pr.array
    ystp_pr = grids3d.arrays.pressure.ystp_pr.array
    zstp_pr = grids3d.arrays.pressure.zstp_pr.array
    spacing_loop!(xnum, gridx_b, xstp_pr)
    spacing_loop!(ynum, gridy_b, ystp_pr)
    spacing_loop!(znum, gridz_b, zstp_pr)
    return nothing
end

function spacing_loop!(
    num::Int64, 
    grid_b::Vector{Float64}, 
    stp_pr::Vector{Float64}
)::Nothing
    for i in 1:(num-2)
        stp_pr[i] = (grid_b[i+2] - grid_b[i])/2.0
    end
    return nothing
end

function calculate_grid_coordinates!(grids::Grids)
    xnum = grids.parameters.geometry.xnum.value
    ynum = grids.parameters.geometry.ynum.value
    gridx_b = grids.arrays.basic.gridx_b.array
    gridy_b = grids.arrays.basic.gridy_b.array
    xstp_pr = grids.arrays.pressure.xstp_pr.array
    ystp_pr = grids.arrays.pressure.ystp_pr.array
    gridx_pr = grids.arrays.pressure.gridx_pr.array
    gridy_pr = grids.arrays.pressure.gridy_pr.array
    coordinate_loop!(xnum, gridx_b, gridx_pr, xstp_pr)
    coordinate_loop!(ynum, gridy_b, gridy_pr, ystp_pr)
end

function calculate_grid_coordinates!(grids3d::Grids3d)
    xnum = grids3d.parameters.geometry.xnum.value
    ynum = grids3d.parameters.geometry.ynum.value
    znum = grids3d.parameters.geometry.znum.value
    xstp_pr = grids3d.arrays.pressure.xstp_pr.array
    ystp_pr = grids3d.arrays.pressure.ystp_pr.array
    zstp_pr = grids3d.arrays.pressure.zstp_pr.array
    gridx_b = grids3d.arrays.basic.gridx_b.array
    gridy_b = grids3d.arrays.basic.gridy_b.array
    gridz_b = grids3d.arrays.basic.gridz_b.array
    gridx_pr = grids3d.arrays.pressure.gridx_pr.array
    gridy_pr = grids3d.arrays.pressure.gridy_pr.array
    gridz_pr = grids3d.arrays.pressure.gridz_pr.array
    coordinate_loop!(xnum, gridx_b, gridx_pr, xstp_pr)
    coordinate_loop!(ynum, gridy_b, gridy_pr, ystp_pr)
    coordinate_loop!(znum, gridz_b, gridz_pr, zstp_pr)
end

function coordinate_loop!(
    num::Int64,
    grid_b::Vector{Float64},
    grid_pr::Vector{Float64},
    stp_pr::Vector{Float64}
)
    grid_pr[1] = (grid_b[1] + grid_b[2])/2.0
    for i in 2:(num-1)
        grid_pr[i] = grid_pr[i-1] + stp_pr[i-1]
    end
end

end # module 