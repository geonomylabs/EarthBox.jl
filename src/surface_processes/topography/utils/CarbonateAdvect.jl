module CarbonateAdvect

import EarthBox.ModelDataContainer: ModelData

""" Vertically advect topography surface based on carbonate growth model.

If iuse_carb is 0, the topography is advected vertically based on erosion and 
sedimentation rates without the carbonate growth model. If iuse_carb is 1, the 
topography is advected vertically based on erosion and sedimentation rates, and 
carbonate growth model.

# Updated Arrays
- `gridt[1,N]`: Elevation of topography node.
"""
function advect_vertically_carbonate!(model::ModelData)::Nothing
    iuse_extrusion = model.melting.parameters.extrusion.iuse_extrusion.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    sec_per_Myr = model.conversion.parameters.sec_per_Myr.value
    time_Myr = timesum/sec_per_Myr

    toponum = model.topography.parameters.topo_grid.toponum.value
    gridt = model.topography.arrays.gridt.array
    erosion_rate = model.topography.parameters.depo_and_erosion_rates.erosion_rate.value
    sedimentation_rate = model.topography.parameters.depo_and_erosion_rates.sedimentation_rate.value
    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value
    grid_carb = model.carbonate.arrays.grid_carb.array
    growth_rad = model.carbonate.parameters.carb_growth_rad.value
    base_rate = model.carbonate.parameters.carb_base_rate.value
    carb_time_myr = model.carbonate.parameters.carb_time_myr.value
    carb_jump_time_myr = model.carbonate.parameters.carb_jump_time_myr.value
    base_rate_jump = model.carbonate.parameters.carb_base_rate_jump.value

    thick_air = model.geometry.parameters.sticky_air_geometry.thick_air.value
    photic_thick_m = model.carbonate.parameters.photic_thick_m.value
    zmax_carb = thick_air + photic_thick_m

    # Setup for carbonate growth
    iuse_carb = 1
    # Advect topography vertically
    # Calculate carb nucleation flags (pvals)
    # Update carb nucleation flags
    if iuse_carb == 1
        if time_Myr < carb_jump_time_myr
            ct_Myr = carb_time_myr
            br = base_rate
        else
            ct_Myr = carb_jump_time_myr
            br = base_rate_jump
        end
        carb_nucleation_flags(model, br, zmax_carb, time_Myr, ct_Myr)
    end

    for j in 1:toponum
        for i in 1:2
            xt = gridt[i, j]
            yt = gridt[i+1, j]
            if iuse_carb == 0
                if yt < y_sealevel
                    v2 = erosion_rate
                else
                    v2 = -sedimentation_rate
                end
            else
                if yt < y_sealevel
                    v2 = erosion_rate
                else
                    if yt < zmax_carb
                        # if surface is in photic zone, find the closest carbonate particle
                        xmin_carb = get_xmin_carb_v2(grid_carb, toponum, xt, zmax_carb)
                        if xmin_carb <= growth_rad
                            pval = 1.0
                        else
                            pval = 0.0
                        end
                        v2 = -sedimentation_rate * pval
                    else
                        v2 = 0.0
                    end
                end
            end
            if iuse_extrusion == 1
                v2 = v2 - gridt[i+5, j]/timestep
            end
            vyt = gridt[i+3, j] + v2
            eo = gridt[i+1, j]
            gridt[i+1, j] = eo + vyt * timestep
        end
    end
    return nothing
end

""" Calculate carbonate nucleation flags.

# Updated Arrays
- `grid_carb[6, N]`: Carbonate array with nucleation probability
"""
function carb_nucleation_flags(
    model::ModelData,
    base_rate::Float64,
    zmax_carb::Float64,
    time_Myr::Float64,
    carb_time_myr::Float64
)::Nothing
    growth_rad = model.carbonate.parameters.carb_growth_rad.value
    grid_carb = model.carbonate.arrays.grid_carb.array
    toponum = model.topography.parameters.topo_grid.toponum.value
    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value

    for j in 1:toponum
        for i in 1:2
            xc = grid_carb[i, j]
            yc = grid_carb[i+1, j]
            if yc >= y_sealevel
                xmin_carb = 1e32
                pp = -9999.0
                pval = -9999.0
                if yc < zmax_carb
                    xmin_carb = get_xmin_carb_v2(grid_carb, toponum, xc, zmax_carb)
                    if xmin_carb <= growth_rad
                        pp = 1.0
                    else
                        pp = base_rate
                    end
                    if time_Myr >= carb_time_myr
                        array1 = zeros(2)
                        array1[1] = 1 - pp
                        array1[2] = pp
                        pval = rand(array1)
                    else
                        pval = 0.0
                    end
                    grid_carb[i+5, j] = pval
                end
            end
        end
    end
    return nothing
end

""" Get distance of closest carbonate grid node in photic zone.

# Returns
- `xmin_carb`: Distance of closest carbonate grid node in photic zone.
"""
function get_xmin_carb_v2(
    grid_carb::Array{Float64,2}, 
    toponum::Int, xt::Float64, 
    zmax_carb::Float64
)::Float64
    xmin_carb = 1e32
    for j in 1:toponum
        for i in 1:2
            xc = grid_carb[i, j]
            yc = grid_carb[i+1, j]
            p = grid_carb[i+5, j]
            # if the particle is carbonate sediment and is in the photic zone
            if yc <= zmax_carb && p > 0.0
                dist = abs(xt - xc)
                if dist < xmin_carb
                    xmin_carb = dist
                end
            end
        end
    end
    return xmin_carb
end

""" Get distance of closest carbonate particle in photic zone.
"""
function get_xmin_carb_using_markers(
    marknum::Int, 
    marker_matid::Array{Int16,1}, 
    marker_y::Array{Float64,1}, 
    marker_x::Array{Float64,1}, 
    xt::Float64, zmax_carb::Float64
)::Float64
    xmin_carb = 1e32
    for mm1 in 1:marknum
        mitype = marker_matid[mm1]
        myy = marker_y[mm1]
        # if the particle is carbonate sediment and is in the photic zone
        if mitype == 2 && myy <= zmax_carb
            mxx = marker_x[mm1]
            dist = abs(xt - mxx)
            if dist < xmin_carb
                xmin_carb = dist
            end
        end
    end
    return xmin_carb
end

""" Advect carbonate surface nodes horizontally.

# Updated Arrays
- `grid_carb[0,N]`: x-location (m) of carbonate node.
- `grid_carb[1,N]`: Elevation (m) of carbonate node.
"""
function advect_carb!(model::ModelData)::Nothing
    toponum = model.topography.parameters.topo_grid.toponum.value
    grid_carb = model.carbonate.arrays.grid_carb.array
    gridt = model.topography.arrays.gridt.array
    dx_topo = model.topography.parameters.topo_grid.dx_topo.value
    timestep = model.timestep.parameters.main_time_loop.timestep.value

    for j in 1:toponum
        for i in 1:2
            vx = grid_carb[i+2, j]
            grid_carb[i, j] = grid_carb[i, j] + vx * timestep
        end
    end

    for j in 1:toponum
        for i in 1:2
            xc = grid_carb[i, j]
            arg = (xc - gridt[i, 1])/dx_topo - 0.5
            xn = floor(Int, arg)
            if xn < 1
                xn = 1
            end
            if xn > toponum - 1
                xn = toponum - 1
            end
            # Compute relative distance to topography node
            dx = (xc - gridt[i, xn])/dx_topo
            # Compute topography elevation above the carbonate node
            dy = gridt[i+1, xn] * (1.0 - dx) + gridt[i+1, xn+1] * dx
            grid_carb[i+1, j] = dy
        end
    end
    return nothing
end

end # module 