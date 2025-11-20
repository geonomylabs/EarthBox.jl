module Dimensions

import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d

function update_grid_dimensions!(grids::Grids)::Nothing
    gridx_b = grids.arrays.basic.gridx_b.array
    gridy_b = grids.arrays.basic.gridy_b.array

    geometry = grids.parameters.geometry
    xnum = geometry.xnum.value
    ynum = geometry.ynum.value

    geometry.xmin.value = gridx_b[1]
    geometry.xmax.value = gridx_b[xnum]
    geometry.ymin.value = gridy_b[1]
    geometry.ymax.value = gridy_b[ynum]
    return nothing
end

function update_grid_dimensions!(grids3d::Grids3d)::Nothing
    gridx_b = grids3d.arrays.basic.gridx_b.array
    gridy_b = grids3d.arrays.basic.gridy_b.array
    gridz_b = grids3d.arrays.basic.gridz_b.array

    geometry = grids3d.parameters.geometry
    xnum = geometry.xnum.value
    ynum = geometry.ynum.value
    znum = geometry.znum.value

    geometry.xmin.value = gridx_b[1]
    geometry.xmax.value = gridx_b[xnum]
    geometry.ymin.value = gridy_b[1]
    geometry.ymax.value = gridy_b[ynum]
    geometry.zmin.value = gridz_b[1]
    geometry.zmax.value = gridz_b[znum]
    return nothing
end

end # module Dimensions
