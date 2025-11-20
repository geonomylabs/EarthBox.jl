module BenchmarkFuncs

import EarthBox.ModelDataContainer: ModelData
import EarthBox.Interpolation.MarkerGridMapping: get_index_of_left_node

function stokes_benchmarks_post_processing(model::ModelData)
    bm_params = model.benchmarks.parameters
    iuse_viscous_block_processing = bm_params.iuse_viscous_block_processing.value

    vel_option = model.stokes_continuity.parameters.velocity_calc_option
    itype_velocity = vel_option.itype_velocity.value

    if itype_velocity == 0
        if iuse_viscous_block_processing == 1
            vis_block_func(model)
        end
    end
end

""" Calculate parameters for the viscous block benchmark.
"""
function vis_block_func(model::ModelData)
    xn, yn, vy_block = calc_vis_block(model)
    cm_yr2m_s = model.conversion.parameters.cm_yr2m_s.value
    println("******************************************************************")
    println("Vy_BLOCK : xn, yn, (m/s), (cm/yr): ", xn, " ", yn, " ", vy_block, " ", 
            vy_block/cm_yr2m_s)
    println("******************************************************************")
end

""" Calculate parameters for the viscous block benchmark.

Returns
-------
- xn : Float64
    x-index of block center
- yn : Float64
    y-index of block center
- vy_block : Float64
    Vertical velocity at block center
"""
function calc_vis_block(model::ModelData)::Tuple{Float64, Float64, Float64}
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    xstpavg = model.grids.parameters.geometry.xstpavg.value
    ystpavg = model.grids.parameters.geometry.ystpavg.value

    gridx_b = model.grids.arrays.basic.gridx_b.array
    gridy_b = model.grids.arrays.basic.gridy_b.array

    vel_arrays = model.stokes_continuity.arrays.staggered_grid_velocity
    vy1 = vel_arrays.vy1.array

    # Define vertical velocity for the initial center of the block
    # Initial position of block center in vy grid
    xblock = (xsize + xstpavg)/2.0
    yblock = 0.2*ysize
    xn = get_index_of_left_node(xblock, gridx_b, xnum)
    yn = get_index_of_left_node(yblock, gridy_b, ynum)
    dx = (xblock - gridx_b[xn])/xstpavg
    dy = (yblock - gridy_b[yn])/ystpavg
    vy_block = (1.0 - dx)*(1.0 - dy)*vy1[yn, xn] +
               dx*(1.0-dy)*vy1[yn, xn+1] +
               (1.0-dx)*dy*vy1[yn+1, xn] +
               dx*dy*vy1[yn+1, xn+1]
    return xn, yn, vy_block
end

""" Call post-processing calculations for convection in a box benchmark.
"""
function conbox_calculations(model::ModelData)
    if model.benchmarks.parameters.iuse_conbox_post_processing.value == 1
        vrms, tavr, nnus1 = conbox_func(model)

        sec_per_myr = model.conversion.parameters.sec_per_Myr.value
        print_conbox_info(
            model.timestep.parameters.main_time_loop.ntimestep.value,
            model.timestep.parameters.main_time_loop.timesum.value/sec_per_myr,
            vrms,
            tavr,
            nnus1
        )
    end
end

""" Perform post-processing calculations for convection in a box benchmark.

Returns
-------
- vrms : Float64
    Root mean square of velocity in m/s
- tavr : Float64
    Average temperature in Kelvins
- nnus1 : Float64
    Nusselt Number
"""
function conbox_func(model::ModelData)
    temperature_bottom = model.bcs.parameters.temperature.temperature_bottom.value
    temperature_top = model.bcs.parameters.temperature.temperature_top.value
    xnum = model.grids.parameters.geometry.xnum.value
    ynum = model.grids.parameters.geometry.ynum.value
    xsize = model.grids.parameters.geometry.xsize.value
    ysize = model.grids.parameters.geometry.ysize.value
    velocity = model.stokes_continuity.arrays.staggered_grid_velocity
    vx1 = velocity.vx1.array
    vy1 = velocity.vy1.array
    tk1 = model.heat_equation.arrays.temperature.tk1.array
    xstp_b = model.grids.arrays.basic.xstp_b.array
    ystp_b = model.grids.arrays.basic.ystp_b.array
    mat_rho = model.materials.arrays.mat_rho.array
    mat_cp = model.materials.arrays.mat_cp.array
    mat_kt = model.materials.arrays.mat_kt.array
    
    vrms = 0.0
    tavr = 0.0
    
    # Note: Changed loop order for column-major Julia arrays
    for j in 1:xnum-1
        for i in 1:ynum-1
            vrms += (vx1[i+1, j]^2 + vx1[i+1, j+1]^2 + 
                    vy1[i, j+1]^2 + vy1[i+1, j+1]^2)/2.0 * 
                    xstp_b[j]*ystp_b[i]
            tavr += (tk1[i, j] + tk1[i+1, j] + 
                    tk1[i, j+1] + tk1[i+1, j+1])/4.0 * 
                    xstp_b[j]*ystp_b[i]
        end
    end
    
    rho = mat_rho[1, 1]
    cp = mat_cp[1]
    k = mat_kt[1, 1]
    vrms = sqrt(vrms/xsize/ysize)*ysize/(k/cp/rho)
    tavr = tavr/xsize/ysize
    nnus1 = 0.0
    
    for j in 1:xnum-1
        dtdy1 = (tk1[1, j] - tk1[2, j] + 
                tk1[1, j+1] - tk1[2, j+1])/2.0/ystp_b[1]
        nnus1 -= dtdy1*xstp_b[j]
    end
    
    nnus1 = ysize*nnus1/((temperature_bottom - temperature_top)*xsize)
    return vrms, tavr, nnus1
end

""" Print information for convection in a box benchmark.
"""
function print_conbox_info(
    ntimestep::Int, 
    timesum::Float64, 
    vrms::Float64, 
    tavr::Float64, nnus1::Float64
)
    println("CONBOX: ntimestep, timesum_Myr, vrms_m_s, tavr_K, nnus1 : ",
            ntimestep, " ", timesum, " ", vrms, " ", tavr, " ", nnus1)
end

end # module 