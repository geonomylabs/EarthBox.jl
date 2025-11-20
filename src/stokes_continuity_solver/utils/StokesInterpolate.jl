module StokesInterpolate

import EarthBox.ModelDataContainer: ModelData

""" Interpolate velocity from staggered grid to basic grid.

# Updated Arrays
- `model.stokes_continuity.arrays.basic_grid_velocity`
  - `vxb.array`: Array((ynum,xnum), Float64)
    x-component of velocity on basic grid.
  - `vyb.array`: Array((ynum,xnum), Float64)
    y-component of velocity on basic grid.

Note: This was located right before the function that advances model time.
"""
function interpolate_stag_velocity_to_basic_grid!(model::ModelData)
    vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array

    vxb = model.stokes_continuity.arrays.basic_grid_velocity.vxb.array
    vyb = model.stokes_continuity.arrays.basic_grid_velocity.vyb.array

    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    # Note: In Julia, we use column-major order, so outer loop is over columns (j)
    for j in 1:xnum
        for i in 1:ynum
            vxb[i, j] = (vx1[i, j] + vx1[i+1, j]) / 2.0
            vyb[i, j] = (vy1[i, j] + vy1[i, j+1]) / 2.0
        end
    end
end

end # module 