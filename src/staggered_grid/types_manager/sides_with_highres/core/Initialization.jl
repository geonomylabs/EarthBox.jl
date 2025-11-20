module Initialization

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import ....Dimensions: update_grid_dimensions!
import ....Spacing: update_average_spacing_of_basic_grid!

""" Initialize basic grid with high resolution along sides.
"""
function initialize!(grids::Grids)
    update_average_spacing_of_basic_grid!(grids)
    calculate_grid_coordinates!(grids)
    update_grid_dimensions!(grids)
end

""" Calculate coordinates for 51x51 basic grid with high-resolution along sides.
"""
function calculate_grid_coordinates!(grids::Grids)
    gridx_b = grids.arrays.basic.gridx_b.array
    gridy_b = grids.arrays.basic.gridy_b.array
    xnum = grids.parameters.geometry.xnum.value
    ynum = grids.parameters.geometry.ynum.value
    xsize = grids.parameters.geometry.xsize.value
    ysize = grids.parameters.geometry.ysize.value

    # Horizontal direction
    for i in 2:11
        gridy_b[i] = gridy_b[i-1] + 10000.0
    end
    for i in 12:16
        gridy_b[i] = gridy_b[i-1] + 10000.0*1.2407234^(convert(Float64, i - 10))
    end
    for i in 17:36
        gridy_b[i] = gridy_b[i-1] + 30000.0
    end
    for i in 37:41
        gridy_b[i] = gridy_b[i-1] + 10000.0*1.2407234^(convert(Float64, 41-i))
    end
    for i in 42:ynum
        gridy_b[i] = gridy_b[i-1] + 10000.0
    end
    gridy_b[ynum] = ysize

    # Vertical direction
    for i in 1:xnum
        gridx_b[i] = gridy_b[i]*xsize/ysize
    end
end

end # module

