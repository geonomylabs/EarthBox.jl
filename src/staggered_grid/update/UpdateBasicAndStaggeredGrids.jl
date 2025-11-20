module UpdateBasicAndStaggeredGrids

import EarthBox.Arrays.ArrayUtils: setzeros!
import EarthBox.ModelDataContainer.Grids2dContainer: Grids
import EarthBox.ModelDataContainer.Grids3dContainer: Grids3d
import ..Spacing: update_spacing_of_basic_grid!
import ..Spacing: calculate_average_grid_spacing_using_grid_arrays
import ..UpdateStaggeredGrids

function update!(grids::Grids)
    setzeros!(grids.arrays.basic.xstp_b)
    setzeros!(grids.arrays.basic.ystp_b)
    update_spacing_of_basic_grid!(grids)

    (
        grids.parameters.geometry.xstpavg.value,
        grids.parameters.geometry.ystpavg.value
    ) = calculate_average_grid_spacing_using_grid_arrays(grids)

    setzeros!(grids.arrays.pressure.xstp_pr)
    setzeros!(grids.arrays.pressure.ystp_pr)
    setzeros!(grids.arrays.staggered_vy.xstp_vy)
    setzeros!(grids.arrays.staggered_vx.ystp_vx)
    UpdateStaggeredGrids.update_staggered_grids!(grids)
end

function update!(grids3d::Grids3d)
    setzeros!(grids3d.arrays.basic.xstp_b)
    setzeros!(grids3d.arrays.basic.ystp_b)
    setzeros!(grids3d.arrays.basic.zstp_b)
    update_spacing_of_basic_grid!(grids3d)

    (
        grids3d.parameters.geometry.xstpavg.value,
        grids3d.parameters.geometry.ystpavg.value,
        grids3d.parameters.geometry.zstpavg.value
    ) = calculate_average_grid_spacing_using_grid_arrays(grids3d)

    setzeros!(grids3d.arrays.pressure.xstp_pr)
    setzeros!(grids3d.arrays.pressure.ystp_pr)
    setzeros!(grids3d.arrays.pressure.zstp_pr)
    
    setzeros!(grids3d.arrays.staggered_vy.xstp_vy)
    setzeros!(grids3d.arrays.staggered_vy.ystp_vy)
    setzeros!(grids3d.arrays.staggered_vy.zstp_vy)

    setzeros!(grids3d.arrays.staggered_vx.xstp_vx)
    setzeros!(grids3d.arrays.staggered_vx.ystp_vx)
    setzeros!(grids3d.arrays.staggered_vx.zstp_vx)

    setzeros!(grids3d.arrays.staggered_vz.xstp_vz)
    setzeros!(grids3d.arrays.staggered_vz.ystp_vz)
    setzeros!(grids3d.arrays.staggered_vz.zstp_vz)

    setzeros!(grids3d.arrays.staggered_sxy.xstp_sxy)
    setzeros!(grids3d.arrays.staggered_sxy.ystp_sxy)
    setzeros!(grids3d.arrays.staggered_sxy.zstp_sxy)

    setzeros!(grids3d.arrays.staggered_sxz.xstp_sxz)
    setzeros!(grids3d.arrays.staggered_sxz.ystp_sxz)
    setzeros!(grids3d.arrays.staggered_sxz.zstp_sxz)

    setzeros!(grids3d.arrays.staggered_syz.xstp_syz)
    setzeros!(grids3d.arrays.staggered_syz.ystp_syz)
    setzeros!(grids3d.arrays.staggered_syz.zstp_syz)

    UpdateStaggeredGrids.update_staggered_grids!(grids3d)
end

end # module 