module TransferBC

import EarthBox.ModelDataContainer: ModelData

""" Transfer boundary conditions to component arrays.

Component arrays are used to build the system of equations.

Updated Arrays
-------------

model.bcs.arrays.vel_comp
    btopx.array: Array((xnum+1,2), Float64)
        Velocity boundary conditions along top boundary:
            vx[1,j] = btopx[j,1] + vx[2,j]*btopx[j,2]

    btopy.array: Array((xnum+1,2), Float64)
        Velocity boundary conditions along top boundary:
            vy[1,j] = btopy[j,1] + vy[2,j]*btopy[j,2]

    bbottomx.array: Array((xnum+1,2), Float64)
        Velocity boundary conditions along bottom boundary:
            vx[ynum+1,j] = bbottomx[j,1] + vx[ynum,j]*bbottomx[j,2]

    bbottomy.array: Array((xnum+1,2), Float64)
        Velocity boundary conditions along bottom boundary:
            vy[ynum+1,j] = bbottomy[j,1] + vy[ynum,j]*bbottomy[j,2]

    bleftx.array: Array((ynum+1,2), Float64)
        Velocity boundary conditions along left boundary:
            vx[i,1] = bleftx[i,1] + vx[i,2]*bleftx[i,2]

    blefty.array: Array((ynum+1,2), Float64)
        Velocity boundary conditions along left boundary:
            vy[i,1] = blefty[i,1] + vy[i,2]*blefty[i,2]

    brightx.array: Array((ynum+1,2), Float64)
        Velocity boundary conditions along right boundary:
            vx[i,xnum+1] = brightx[i,1] + vx[i,xnum]*brightx[i,2]

    brighty.array: Array((ynum+1,2), Float64)
        Velocity boundary conditions along right boundary:
            vy[i,xnum+1] = brighty[i,1] + vx[i,xnum]*brighty[i,2]
"""
function transfer_velocity_bcs_to_component_arrays!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    btop = model.bcs.arrays.velocity.btop.array
    bbottom = model.bcs.arrays.velocity.bbottom.array
    bleft = model.bcs.arrays.velocity.bleft.array
    bright = model.bcs.arrays.velocity.bright.array

    btopx = model.bcs.arrays.vel_comp.btopx.array
    bbottomx = model.bcs.arrays.vel_comp.bbottomx.array
    bleftx = model.bcs.arrays.vel_comp.bleftx.array
    brightx = model.bcs.arrays.vel_comp.brightx.array

    btopy = model.bcs.arrays.vel_comp.btopy.array
    bbottomy = model.bcs.arrays.vel_comp.bbottomy.array
    blefty = model.bcs.arrays.vel_comp.blefty.array
    brighty = model.bcs.arrays.vel_comp.brighty.array

    for i in 1:xnum+1
        btopx[i,1] = btop[i,1]
        btopx[i,2] = btop[i,2]
        btopy[i,1] = btop[i,3]
        btopy[i,2] = btop[i,4]
        bbottomx[i,1] = bbottom[i,1]
        bbottomx[i,2] = bbottom[i,2]
        bbottomy[i,1] = bbottom[i,3]
        bbottomy[i,2] = bbottom[i,4]
    end

    for i in 1:ynum+1
        bleftx[i,1] = bleft[i,1]
        bleftx[i,2] = bleft[i,2]
        blefty[i,1] = bleft[i,3]
        blefty[i,2] = bleft[i,4]
        brightx[i,1] = bright[i,1]
        brightx[i,2] = bright[i,2]
        brighty[i,1] = bright[i,3]
        brighty[i,2] = bright[i,4]
    end

    return nothing
end

end # module