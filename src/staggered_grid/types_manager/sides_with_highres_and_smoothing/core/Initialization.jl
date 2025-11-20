module Initialization

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import ....Dimensions: update_grid_dimensions!
import ....Spacing: update_average_spacing_of_basic_grid!

""" Initialize grid with high resolution at sides and smoothing.
"""
function initialize!(grids::Grids)
    update_average_spacing_of_basic_grid!(grids)
    calculate_grid_coordinates!(grids)
    update_grid_dimensions!(grids)
end

function calculate_grid_coordinates!(grids::Grids)
    gridx_b = grids.arrays.basic.gridx_b.array
    gridy_b = grids.arrays.basic.gridy_b.array
    xnum = grids.parameters.geometry.xnum.value
    ynum = grids.parameters.geometry.ynum.value
    xsize = grids.parameters.geometry.xsize.value
    ysize = grids.parameters.geometry.ysize.value

    # Number of nodes associated with thermal boundary layer
    ynum_BL = Int64(ynum/4.0)
    # Number of nodes associated with smoothing zone
    ynum_SM = Int64(ynum_BL/2.0)
    # Thickness of high-resolution zones
    L_BL = ysize/10.0
    # Thickness of intermediate resolution zone
    L_SM = L_BL
    ynum_M = ynum - 2.0*ynum_BL - 2.0*ynum_SM + 4.0

    ip2 = ynum_BL - 1
    ip3 = ip2 + ynum_SM - 1
    ip4 = ip3 + ynum_M - 1
    ip5 = ip4 + ynum_SM - 1

    L_M = ysize - L_BL*2.0 - L_SM*2.0

    ystp_BL = L_BL/convert(Float64, ynum_BL - 1)
    ystp_SM = L_SM/convert(Float64, ynum_SM - 1)
    ystp_M = L_M/convert(Float64, ynum_M - 1)

    # Horizontal direction
    for i in 2:ip2
        gridy_b[i] = gridy_b[i-1] + ystp_BL
    end
    for i in (ip2+1):ip3
        gridy_b[i] = gridy_b[i-1] + ystp_SM
    end
    for i in (ip3+1):ip4
        gridy_b[i] = gridy_b[i-1] + ystp_M
    end
    for i in (ip4+1):ip5
        gridy_b[i] = gridy_b[i-1] + ystp_SM
    end
    for i in (ip5+1):ynum
        gridy_b[i] = gridy_b[i-1] + ystp_BL
    end

    gridy_b[ynum] = ysize

    # Vertical direction
    for i in 1:xnum
        gridx_b[i] = gridy_b[i]*xsize/ysize
    end
end

end # module 