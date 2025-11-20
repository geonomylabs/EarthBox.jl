module Kinematics

import EarthBox.ModelDataContainer: ModelData

"""
    solid_body_rotation!(model::ModelData)

Compute velocity using solid body rotation.

# Updated Arrays
- `vx1`: x-component of velocity defined on staggered vx grid
- `vy1`: y-component of velocity defined on staggered vy grid
"""
function solid_body_rotation!(model::ModelData)::Nothing
    # Get parameters from model
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    xstpavg = model.grids.parameters.geometry.xstpavg.value
    ystpavg = model.grids.parameters.geometry.ystpavg.value
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    vrotation = model.bcs.parameters.velocity.velocity_rotation.value
    vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array

    for j in 1:(xnum+1)
        for i in 1:(ynum+1)
            if j <= xnum
                # Relative distance of vx node from the model center
                dx = (Float64(j-1)*xstpavg - xsize/2.0)/(xsize/2.0)
                dy = ((Float64(i-1)-0.5)*ystpavg - ysize/2.0)/(xsize/2.0)
                vx1[i, j] = -vrotation*dy
            end
            if i <= ynum
                # Relative distance of vy node from the model center
                dx = (Float64(j-1-0.5)*xstpavg - xsize/2.0)/(xsize/2.0)
                dy = (Float64(i-1)*ystpavg - ysize/2.0)/(xsize/2.0)
                vy1[i, j] = vrotation*dx
            end
        end
    end
    return nothing
end

end # module 