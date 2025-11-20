module Initialization

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import ....Dimensions: update_grid_dimensions!
import ....Spacing: update_average_spacing_of_basic_grid!

""" Initialize basic 2D grid with uniform spacing.
"""
function initialize!(grids::Grids)
    update_average_spacing_of_basic_grid!(grids)
    calculate_grid_coordinates!(grids)
    update_grid_dimensions!(grids)
end

""" Initialize basic 3D grid with uniform spacing.
"""
function initialize!(grids3d::Grids3d)
    update_average_spacing_of_basic_grid!(grids3d)
    calculate_grid_coordinates!(grids3d)
    update_grid_dimensions!(grids3d)
end


""" Initialize uniform basic 2D grid.
"""
function calculate_grid_coordinates!(grids::Grids)
    gridx_b = grids.arrays.basic.gridx_b.array
    gridy_b = grids.arrays.basic.gridy_b.array
    xnum = grids.parameters.geometry.xnum.value
    ynum = grids.parameters.geometry.ynum.value
    xstpavg = grids.parameters.geometry.xstpavg.value
    ystpavg = grids.parameters.geometry.ystpavg.value

    gridx_b[1] = 0.0
    for i in 2:xnum
        gridx_b[i] = gridx_b[i-1] + xstpavg
    end
    gridy_b[1] = 0.0
    for i in 2:ynum
        gridy_b[i] = gridy_b[i-1] + ystpavg
    end
end

""" Initialize uniform basic 3D grid.
"""
function calculate_grid_coordinates!(grids3d::Grids3d)
    gridx_b = grids3d.arrays.basic.gridx_b.array
    gridy_b = grids3d.arrays.basic.gridy_b.array
    gridz_b = grids3d.arrays.basic.gridz_b.array
    xnum = grids3d.parameters.geometry.xnum.value
    ynum = grids3d.parameters.geometry.ynum.value
    znum = grids3d.parameters.geometry.znum.value
    xsize = grids3d.parameters.geometry.xsize.value
    ysize = grids3d.parameters.geometry.ysize.value
    zsize = grids3d.parameters.geometry.zsize.value

    xstp = xsize/(xnum - 1)
    ystp = ysize/(ynum - 1)
    zstp = zsize/(znum - 1)

    gridx_b[1] = 0.0
    for i in 2:xnum
        gridx_b[i] = gridx_b[i-1] + xstp
    end
    gridy_b[1] = 0.0
    for i in 2:ynum
        gridy_b[i] = gridy_b[i-1] + ystp
    end
    gridz_b[1] = 0.0
    for i in 2:znum
        gridz_b[i] = gridz_b[i-1] + zstp
    end
end

end # module 