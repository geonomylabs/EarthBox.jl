module GridStrainRate

import EarthBox.ModelDataContainer: ModelData
import EarthBox.PrintFuncs: @timeit_memit

function update_grid_strain_rate_and_spin!(model::ModelData)::Nothing
    @timeit_memit "Finished updating grid strain rate and spin" begin
        calculate_normal_strain_rate!(model)
        calculate_shear_strain_rate!(model)
        calculate_second_invariant!(model)
        calculate_spin!(model)
    end
    return nothing
end

""" Compute normal strain rate on pressure grid.

Updated Arrays
==============
model.stokes_continuity.arrays.strain_rate_and_spin
---------------------------------------------------
exx: Array{Float64,2}
    Normal strain rate (1/s) on pressure grid.

Mathematical Background
======================
The following general equations are used to compute strain rate:

    exx = -eyy = 1/2(dvx/dx - dvy/dy)

where vx and vy are the x- and y-components of velocity, respectively.
"""
function calculate_normal_strain_rate!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    xstp_b = model.grids.arrays.basic.xstp_b.array
    ystp_b = model.grids.arrays.basic.ystp_b.array

    vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array

    model.stokes_continuity.arrays.strain_rate_and_spin.exx.array .= 0.0
    exx = model.stokes_continuity.arrays.strain_rate_and_spin.exx.array

    for j in 1:xnum
        for i in 1:ynum
            if i < ynum && j < xnum
                exx[i,j] = 0.5 * (
                    (vx1[i+1,j+1] - vx1[i+1,j])/xstp_b[j] -
                    (vy1[i+1,j+1] - vy1[i,j+1])/ystp_b[i]
                )
            end
        end
    end
    return nothing
end

""" Calculate shear strain rate on basic grid.

Updated Arrays
==============
model.stokes_continuity.arrays.strain_rate_and_spin
---------------------------------------------------
exy: Array{Float64,2}
    Shear strain rate (1/s) on basic grid.

Mathematical Background
======================
The following general equations are used to compute strain rate:

    exy = eyx = 1/2(dvx/dy + dvy/dx)

where vx and vy are the x- and y-components of velocity, respectively.
"""
function calculate_shear_strain_rate!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array
    ystp_vx = model.grids.arrays.staggered_vx.ystp_vx.array

    vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array

    model.stokes_continuity.arrays.strain_rate_and_spin.exy.array .= 0.0
    exy = model.stokes_continuity.arrays.strain_rate_and_spin.exy.array

    for j in 1:xnum
        for i in 1:ynum
            exy[i,j] = 0.5 * (
                (vx1[i+1,j] - vx1[i,j])/ystp_vx[i] +
                (vy1[i,j+1] - vy1[i,j])/xstp_vy[j]
            )
        end
    end
    return nothing
end

""" Compute second invariant of deviatoric strain rate.

Updated Arrays
==============
model.stokes_continuity.arrays.strain_rate_and_spin
---------------------------------------------------
eii: Array{Float64,2}
    Second invariant of deviatoric strain rate (1/s) on pressure grid.

Mathematical Background
----------------------
The second invariant of the strain rate tensor is calculated using the
following equation:

    eii = (exx^2 + exy^2)^0.5

Note that exx^2 = eyy^2. Therefore, 0.5exx^2 + 0.5eyy^2 = exx^2.
Also, since exy^2 = eyx^2, 0.5exy^2 + 0.5eyx^2 = exy^2.
"""
function calculate_second_invariant!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value

    exx = model.stokes_continuity.arrays.strain_rate_and_spin.exx.array
    exy = model.stokes_continuity.arrays.strain_rate_and_spin.exy.array

    model.stokes_continuity.arrays.strain_rate_and_spin.eii.array .= 0.0
    eii = model.stokes_continuity.arrays.strain_rate_and_spin.eii.array

    for j in 1:xnum
        for i in 1:ynum
            if i > 1 && j > 1
                eii[i-1,j-1] = sqrt(
                    exx[i-1,j-1]^2.0 +
                    (exy[i-1,j-1]^2.0 + exy[i,j-1]^2.0 +
                     exy[i-1,j]^2.0 + exy[i,j]^2.0)/4.0
                )
            end
        end
    end
    return nothing
end

""" Compute spin.

Updated Arrays
==============
model.stokes_continuity.arrays.strain_rate_and_spin
---------------------------------------------------
esp: Array{Float64,2}
    Spin (1/s) on basic grid.

Mathematical Background
======================
The following equation is used to compute spin:

    esp = 1/2(dvy/dx - dvx/dy)

where vx and vy are the x- and y-components of velocity, respectively.

Spin is positive for clockwise rotation when x-axis is directed rightward
and y-axis is directed downward.
"""
function calculate_spin!(model::ModelData)::Nothing
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    xstp_vy = model.grids.arrays.staggered_vy.xstp_vy.array
    ystp_vx = model.grids.arrays.staggered_vx.ystp_vx.array

    vx1 = model.stokes_continuity.arrays.staggered_grid_velocity.vx1.array
    vy1 = model.stokes_continuity.arrays.staggered_grid_velocity.vy1.array

    model.stokes_continuity.arrays.strain_rate_and_spin.esp.array .= 0.0
    esp = model.stokes_continuity.arrays.strain_rate_and_spin.esp.array

    for j in 1:xnum
        for i in 1:ynum
            esp[i,j] = 0.5 * (
                (vy1[i,j+1] - vy1[i,j])/xstp_vy[j] -
                (vx1[i+1,j] - vx1[i,j])/ystp_vx[i]
            )
        end
    end
    return nothing
end

end # module 