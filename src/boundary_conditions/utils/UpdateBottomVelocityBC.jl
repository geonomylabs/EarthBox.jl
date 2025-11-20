module UpdateBottomVelocityBC

import EarthBox.ModelDataContainer: ModelData

""" Update bottom boundary condition for extension model with moving boundaries.

New model width and height requires that the y-component of velocity
is re-calculated. Therefore, the bottom velocity boundary condition needs
to be re-calculated. The boundary condition is free slip with prescribed
inward velocity.

# Updated Arrays
model.bcs.arrays.velocity
    `bbottom.array::Matrix{Float64}` (xnum+1,4):
        Velocity boundary conditions along bottom boundary:
            vx[ynum+1,j] = bbottom[j,1] + vx[ynum,j]*bbottom[j,2]
            vy[ynum+1,j] = bbottom[j,3] + vy[ynum,j]*bbottom[j,4]

model.bcs.arrays.vel_comp
    `bbottomy.array::Matrix{Float64}` (xnum+1,2):
        Velocity boundary conditions along bottom boundary:
            vy[ynum+1,j] = bbottomy[j,1] + vy[ynum,j]*bbottomy[j,2]
"""
function update_bc(model::ModelData)::Nothing
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    xnum = model.grids.parameters.geometry.xnum.value
    full_velocity_extension = model.bcs.parameters.velocity.full_velocity_extension.value

    bbottom = model.bcs.arrays.velocity.bbottom.array
    bbottomy = model.bcs.arrays.vel_comp.bbottomy.array

    for j in 1:xnum+1
        bbottom[j,3] = -full_velocity_extension/xsize*ysize
        bbottom[j,4] = 0.0
        bbottomy[j,1] = -full_velocity_extension/xsize*ysize
        bbottomy[j,2] = 0.0
    end

    return nothing
end

end # module UpdateBottomVelocityBC 