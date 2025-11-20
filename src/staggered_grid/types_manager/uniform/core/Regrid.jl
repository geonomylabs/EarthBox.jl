module Regrid

import EarthBox.ModelDataContainer: ModelData

"""
    regrid_uniform_grid_using_extension_velocity!(model::ModelData)

Regrid uniform grid using the extension velocity.

Update Arrays
------------
grids.arrays.basic
    gridx_b.array: Array((xnum), Float64)
        X-coordinates of basic grid in meters.

    gridy_b.array: Array((ynum), Float64)
        Y-coordinates of basic grid in meters.
"""
function regrid_uniform_grid_using_extension_velocity!(
    model::ModelData
)::Nothing
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    full_velocity_extension = model.bcs.parameters.velocity.full_velocity_extension.value
    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    xstpavg = model.grids.parameters.geometry.xstpavg.value
    ystpavg = model.grids.parameters.geometry.ystpavg.value

    gridx_b[1] = gridx_b[1] - full_velocity_extension/2.0*timestep
    for i in 2:xnum
        gridx_b[i] = gridx_b[i-1] + xstpavg
    end
    gridy_b[1] = 0.0
    for i in 2:ynum
        gridy_b[i] = gridy_b[i-1] + ystpavg
    end
    return nothing
end

end # module
