module Residuals

import ...LevelManager: LevelData, LevelData2d
import ...Domain: on_vx_boundary3d, on_vy_boundary3d, on_vz_boundary3d
import ...Domain: on_vx_boundary2d, on_vy_boundary2d
import ...ArrayStats

function compute_residuals!(
    level_data::LevelData
)::Tuple{Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}, Array{Float64, 3}}
    vx = level_data.vx.array
    vy = level_data.vy.array
    vz = level_data.vz.array
    pr = level_data.pr.array

    xnum = level_data.grid.parameters.geometry.xnum.value
    ynum = level_data.grid.parameters.geometry.ynum.value
    znum = level_data.grid.parameters.geometry.znum.value

    # Initialize residual arrays
    ΔRx = zeros(Float64, size(vx))
    ΔRy = zeros(Float64, size(vy))
    ΔRz = zeros(Float64, size(vz))
    ΔRc = zeros(Float64, size(pr))

    # Computing final state of residuals
    Threads.@threads for k = 1:znum+1
        for j = 1:xnum+1
            for i = 1:ynum+1
                if j < xnum+1
                    if on_vx_boundary3d(i, j, k, ynum, xnum, znum)
                        ΔRx[i,j,k] = 0.0
                    else
                        ΔRx[i,j,k], _ = calculate_vx_residual(i, j, k, level_data)
                    end
                end
                if i < ynum+1
                    if on_vy_boundary3d(i, j, k, ynum, xnum, znum)
                        ΔRy[i,j,k] = 0.0
                    else
                        ΔRy[i,j,k], _ = calculate_vy_residual(i, j, k, level_data)
                    end
                end
                if k < znum+1
                    if on_vz_boundary3d(i, j, k, ynum, xnum, znum)
                        ΔRz[i,j,k] = 0.0
                    else
                        ΔRz[i,j,k], _ = calculate_vz_residual(i, j, k, level_data)
                    end
                end
                if i < ynum && j < xnum && k < znum
                    ΔRc[i,j,k] = calculate_pressure_residual(i, j, k, level_data)
                end
            end
        end
    end
    return ΔRx, ΔRy, ΔRz, ΔRc
end

function compute_residuals!(
    level_data::LevelData2d
)::Tuple{Array{Float64, 2}, Array{Float64, 2}, Array{Float64, 2}}
    vx = level_data.vx.array
    vy = level_data.vy.array
    pr = level_data.pr.array

    xnum = level_data.grid.parameters.geometry.xnum.value
    ynum = level_data.grid.parameters.geometry.ynum.value

    # Initialize residual arrays
    ΔRx = zeros(Float64, size(vx))
    ΔRy = zeros(Float64, size(vy))
    ΔRc = zeros(Float64, size(pr))

    # Computing final state of residuals
    for j = 1:xnum
        for i = 1:ynum
            if j < xnum
                if on_vx_boundary2d(i, j, ynum, xnum)
                    ΔRx[i,j] = 0.0
                else
                    ΔRx[i,j], _ = calculate_vx_residual(i, j, level_data)
                end
            end
            if i < ynum+1
                if on_vy_boundary2d(i, j, ynum, xnum)
                    ΔRy[i,j] = 0.0
                else
                    ΔRy[i,j], _ = calculate_vy_residual(i, j, level_data)
                end
            end
            if i < ynum && j < xnum
                ΔRc[i,j] = calculate_pressure_residual(i, j, level_data)
            end
        end
    end
    return ΔRx, ΔRy, ΔRc
end

function calculate_vx_residual(
    i::Int64,
    j::Int64,
    k::Int64,
    level_data::LevelData
)::Tuple{Float64, Float64}

    vx = level_data.vx.array
    vy = level_data.vy.array
    vz = level_data.vz.array
    pr = level_data.pr.array
    etaxy = level_data.etaxy.array
    etaxz = level_data.etaxz.array
    etan = level_data.etan.array
    RX = level_data.RX.array

    ystp_b = level_data.grid.arrays.basic.ystp_b.array
    xstp_b = level_data.grid.arrays.basic.xstp_b.array
    zstp_b = level_data.grid.arrays.basic.zstp_b.array
    ystp_vx = level_data.grid.arrays.staggered_vx.ystp_vx.array
    zstp_vx = level_data.grid.arrays.staggered_vx.zstp_vx.array
    xstp_vy = level_data.grid.arrays.staggered_vy.xstp_vy.array

    ΔxL =  xstp_b[j-1]
    ΔxR =  xstp_b[j  ]
    ΔxC = xstp_vy[j  ]

    ΔyU = ystp_vx[i-1]
    ΔyD = ystp_vx[i  ]
    ΔyC =  ystp_b[i-1]

    ΔzF = zstp_vx[k  ]
    ΔzB = zstp_vx[k-1]
    ΔzC =  zstp_b[k-1]

    # Define viscosity, pressure and velocity terms relative to
    # VxC = vx[i,j,k]
    η_nR  = etan[i-1, j  , k-1]
    η_nL  = etan[i-1, j-1, k-1]
    η_xyU = etaxy[i-1, j  , k-1]
    η_xyD = etaxy[i  , j  , k-1]
    η_xzF = etaxz[i-1, j  , k  ]
    η_xzB = etaxz[i-1, j  , k-1]
    pL    = pr[i-1, j-1, k-1]
    pR    = pr[i-1, j  , k-1]
    vxL   = vx[i  , j-1, k  ]
    vxC   = vx[i  , j  , k  ]
    vxR   = vx[i  , j+1, k  ]
    vxU   = vx[i-1, j  , k  ]
    vxD   = vx[i+1, j  , k  ]
    vxF   = vx[i  , j  , k+1]
    vxB   = vx[i  , j  , k-1]
    vyUL  = vy[i-1, j  , k  ]
    vyUR  = vy[i-1, j+1, k  ]
    vyLL  = vy[i  , j  , k  ]
    vyLR  = vy[i  , j+1, k  ]
    vzRF  = vz[i  , j+1, k  ]
    vzLF  = vz[i  , j  , k  ]
    vzRB  = vz[i  , j+1, k-1] 
    vzLB  = vz[i  , j  , k-1]

    dp_dx = (pR - pL)/ΔxC
    dσxx_dx = 2.0*(η_nR*(vxR - vxC)/ΔxR - η_nL*(vxC - vxL)/ΔxL)/ΔxC
    dσxy_dy = (
        η_xyD*((vxD - vxC)/ΔyD + (vyLR - vyLL)/ΔxC) 
        - η_xyU*((vxC - vxU)/ΔyU + (vyUR - vyUL)/ΔxC)
        )/ΔyC
    dσxz_dz = (
        η_xzF*((vzRF - vzLF)/ΔxC + (vxF - vxC)/ΔzF)
        - η_xzB*((vzRB - vzLB)/ΔxC + (vxC - vxB)/ΔzB)
        )/ΔzC

    ΔR = RX[i,j,k] + dp_dx - dσxx_dx - dσxy_dy - dσxz_dz

    # Coefficients for VxC
    Coef_vxC = (
        - 2.0*(η_nR/ΔxR/ΔxC  + η_nL/ΔxL/ΔxC)
        - (η_xyD/ΔyD/ΔyC + η_xyU/ΔyU/ΔyC)
        - (η_xzF/ΔzF/ΔzC + η_xzB/ΔzB/ΔzC)
        )

    return ΔR, Coef_vxC
end

function calculate_vx_residual(
    i::Int64,
    j::Int64,
    level_data::LevelData2d
)::Tuple{Float64, Float64}

    vx = level_data.vx.array
    vy = level_data.vy.array
    pr = level_data.pr.array
    etas = level_data.etas.array
    etan = level_data.etan.array
    RX = level_data.RX.array

    ystp_b = level_data.grid.arrays.basic.ystp_b.array
    xstp_b = level_data.grid.arrays.basic.xstp_b.array
    ystp_vx = level_data.grid.arrays.staggered_vx.ystp_vx.array
    xstp_vy = level_data.grid.arrays.staggered_vy.xstp_vy.array

    ΔxL =  xstp_b[j-1]
    ΔxR =  xstp_b[j  ]
    ΔxC = xstp_vy[j  ]

    ΔyU = ystp_vx[i-1]
    ΔyD = ystp_vx[i  ]
    ΔyC =  ystp_b[i-1]

    # Define viscosity, pressure and velocity terms relative to
    # VxC = vx[i,j]
    η_nR  = etan[i-1, j  ]
    η_nL  = etan[i-1, j-1]
    η_xyU = etas[i-1, j  ]
    η_xyD = etas[i  , j  ]
    pL    =   pr[i-1, j-1]
    pR    =   pr[i-1, j  ]
    vxL   =   vx[i  , j-1]
    vxC   =   vx[i  , j  ]
    vxR   =   vx[i  , j+1]
    vxU   =   vx[i-1, j  ]
    vxD   =   vx[i+1, j  ]
    vyUL  =   vy[i-1, j  ]
    vyUR  =   vy[i-1, j+1]
    vyLL  =   vy[i  , j  ]
    vyLR  =   vy[i  , j+1]

    dp_dx = (pR - pL)/ΔxC
    dσxx_dx = 2.0*(η_nR*(vxR - vxC)/ΔxR - η_nL*(vxC - vxL)/ΔxL)/ΔxC
    dσxy_dy = (
          η_xyD*((vxD - vxC)/ΔyD + (vyLR - vyLL)/ΔxC) 
        - η_xyU*((vxC - vxU)/ΔyU + (vyUR - vyUL)/ΔxC)
        )/ΔyC
    ΔR = RX[i,j] + dp_dx - dσxx_dx - dσxy_dy

    # Coefficients for VxC
    Coef_vxC = (
        - 2.0*(η_nR/ΔxR/ΔxC  + η_nL/ΔxL/ΔxC)
        - (η_xyD/ΔyD/ΔyC + η_xyU/ΔyU/ΔyC)
        )

    return ΔR, Coef_vxC
end

# y-Stokes equation dSIGMAyx/dx+dSIGMAyy/dy+SIGMAyz/dz-dP/dy=RY
function calculate_vy_residual(
    i::Int64,
    j::Int64,
    k::Int64,
    level_data::LevelData
)::Tuple{Float64, Float64}

    vx = level_data.vx.array
    vy = level_data.vy.array
    vz = level_data.vz.array
    pr = level_data.pr.array
    etaxy = level_data.etaxy.array
    etayz = level_data.etayz.array
    etan = level_data.etan.array
    RY = level_data.RY.array

    xstp_b = level_data.grid.arrays.basic.xstp_b.array
    ystp_b = level_data.grid.arrays.basic.ystp_b.array
    zstp_b = level_data.grid.arrays.basic.zstp_b.array
    ystp_vx = level_data.grid.arrays.staggered_vx.ystp_vx.array
    zstp_vx = level_data.grid.arrays.staggered_vx.zstp_vx.array
    xstp_vy = level_data.grid.arrays.staggered_vy.xstp_vy.array

    ΔxL =  xstp_b[j-1]
    ΔxR = xstp_vy[j  ]
    ΔxC = xstp_vy[j-1]

    ΔyU =  ystp_b[i-1]
    ΔyD =  ystp_b[i  ]
    ΔyC = ystp_vx[i  ]

    ΔzF = zstp_vx[k  ]
    ΔzB = zstp_vx[k-1]
    ΔzC =  zstp_b[k-1]

    pD = pr[i  , j-1, k-1]
    pU = pr[i-1, j-1, k-1]

    η_nD  =  etan[i  , j-1, k-1]
    η_nU  =  etan[i-1, j-1, k-1]
    η_xyR = etaxy[i  , j  , k-1]
    η_xyL = etaxy[i  , j-1, k-1]
    η_yzF = etayz[i  , j-1, k  ]
    η_yzB = etayz[i  , j-1, k-1]

    vyU  = vy[i-1, j  , k  ]
    vyC  = vy[i  , j  , k  ]
    vyD  = vy[i+1, j  , k  ]
    vyR  = vy[i  , j+1, k  ]
    vyL  = vy[i  , j-1, k  ]
    vyF  = vy[i  , j  , k+1]
    vyB  = vy[i  , j  , k-1]
    vxLR = vx[i+1, j  , k  ]
    vxUR = vx[i  , j  , k  ]
    vxLL = vx[i+1, j-1, k  ]
    vxUL = vx[i  , j-1, k  ]
    vzDF = vz[i+1, j  , k  ]
    vzUF = vz[i  , j  , k  ]
    vzDB = vz[i+1, j  , k-1]
    vzUB = vz[i  , j  , k-1]

    dp_dy = (pD - pU)/ΔyC
    dσyy_dy = 2.0*(η_nD*(vyD - vyC)/ΔyD - η_nU*(vyC - vyU)/ΔyU)/ΔyC
    dσyx_dx = (
          η_xyR*((vyR - vyC)/ΔxR + (vxLR - vxUR)/ΔyC) 
        - η_xyL*((vyC - vyL)/ΔxL + (vxLL - vxUL)/ΔyC)
        )/ΔxC
    dσyz_dz = (
          η_yzF*((vyF - vyC)/ΔzF + (vzDF - vzUF)/ΔyC)
        - η_yzB*((vyC - vyB)/ΔzB + (vzDB - vzUB)/ΔyC)
        )/ΔzC

    ΔR = RY[i,j,k] + dp_dy - dσyy_dy - dσyx_dx - dσyz_dz

    Coef_vyC = (
        - 2.0*(η_nD/ΔyD/ΔyC +  η_nU/ΔyU/ΔyC)
        -     (η_xyR/ΔxR/ΔxC + η_xyL/ΔxL/ΔxC) 
        -     (η_yzF/ΔzF/ΔzC + η_yzB/ΔzB/ΔzC)
        )

    return ΔR, Coef_vyC
end

function calculate_vy_residual(
    i::Int64,
    j::Int64,
    level_data::LevelData2d
)::Tuple{Float64, Float64}

    vx = level_data.vx.array
    vy = level_data.vy.array
    pr = level_data.pr.array
    etas = level_data.etas.array
    etan = level_data.etan.array
    RY = level_data.RY.array

    xstp_b = level_data.grid.arrays.basic.xstp_b.array
    ystp_b = level_data.grid.arrays.basic.ystp_b.array
    ystp_vx = level_data.grid.arrays.staggered_vx.ystp_vx.array
    xstp_vy = level_data.grid.arrays.staggered_vy.xstp_vy.array

    ΔxL =  xstp_b[j-1]
    ΔxR = xstp_vy[j  ]
    ΔxC = xstp_vy[j-1]

    ΔyU =  ystp_b[i-1]
    ΔyD =  ystp_b[i  ]
    ΔyC = ystp_vx[i  ]

    pD = pr[i  , j-1]
    pU = pr[i-1, j-1]

    η_nD  = etan[i  , j-1]
    η_nU  = etan[i-1, j-1]
    η_xyR = etas[i  , j  ]
    η_xyL = etas[i  , j-1]

    vyU  = vy[i-1, j  ]
    vyC  = vy[i  , j  ]
    vyD  = vy[i+1, j  ]
    vyR  = vy[i  , j+1]
    vyL  = vy[i  , j-1]
    vxLR = vx[i+1, j  ]
    vxUR = vx[i  , j  ]
    vxLL = vx[i+1, j-1]
    vxUL = vx[i  , j-1]

    dp_dy = (pD - pU)/ΔyC
    dσyy_dy = 2.0*(η_nD*(vyD - vyC)/ΔyD - η_nU*(vyC - vyU)/ΔyU)/ΔyC
    dσyx_dx = (
          η_xyR*((vyR - vyC)/ΔxR + (vxLR - vxUR)/ΔyC) 
        - η_xyL*((vyC - vyL)/ΔxL + (vxLL - vxUL)/ΔyC)
        )/ΔxC

    ΔR = RY[i,j] + dp_dy - dσyy_dy - dσyx_dx

    Coef_vyC = (
        - 2.0*(η_nD/ΔyD/ΔyC +  η_nU/ΔyU/ΔyC)
        -     (η_xyR/ΔxR/ΔxC + η_xyL/ΔxL/ΔxC) 
        )

    return ΔR, Coef_vyC
end

# z-Stokes equation dSIGMAzx/dx+dSIGMAzy/dy+SIGMAzz/dz-dP/dz=RZ
function calculate_vz_residual(
    i::Int64,
    j::Int64,
    k::Int64,
    level_data::LevelData
)::Tuple{Float64, Float64}

    vx = level_data.vx.array
    vy = level_data.vy.array
    vz = level_data.vz.array
    pr = level_data.pr.array
    etaxz = level_data.etaxz.array
    etayz = level_data.etayz.array
    etan = level_data.etan.array
    RZ = level_data.RZ.array

    ystp_b = level_data.grid.arrays.basic.ystp_b.array
    xstp_b = level_data.grid.arrays.basic.xstp_b.array
    zstp_b = level_data.grid.arrays.basic.zstp_b.array
    ystp_vx = level_data.grid.arrays.staggered_vx.ystp_vx.array
    zstp_vx = level_data.grid.arrays.staggered_vx.zstp_vx.array
    xstp_vy = level_data.grid.arrays.staggered_vy.xstp_vy.array

    ΔzR =  zstp_b[k  ]
    ΔzL =  zstp_b[k-1]
    ΔzC = zstp_vx[k  ]

    ΔxF = xstp_vy[j  ]
    ΔxB = xstp_vy[j-1]
    ΔxC =  xstp_b[j-1]

    ΔyU = ystp_vx[i-1]
    ΔyD = ystp_vx[i  ]
    ΔyC =  ystp_b[i-1]

    pF = pr[i-1,j-1,k]
    pB = pr[i-1,j-1,k-1]

    η_nF  =  etan[i-1, j-1, k  ]
    η_nB  =  etan[i-1, j-1, k-1]
    η_xzR = etaxz[i-1, j  , k  ]
    η_xzL = etaxz[i-1, j-1, k  ]
    η_yzU = etayz[i-1, j-1, k  ]
    η_yzD = etayz[i  , j-1, k  ]

    vxRF = vx[i  , j  , k+1]
    vxLF = vx[i  , j  , k  ]
    vxRB = vx[i  , j-1, k+1]
    vxLB = vx[i  , j-1, k  ]

    vyLR = vy[i  , j, k+1]
    vyLL = vy[i  , j, k  ]
    vyUR = vy[i-1, j, k+1]
    vyUL = vy[i-1, j, k  ]

    vzL = vz[i  , j  , k-1]
    vzC = vz[i  , j  , k  ]
    vzR = vz[i  , j  , k+1]
    vzF = vz[i  , j+1  , k]
    vzB = vz[i  , j-1  , k]
    vzD = vz[i+1, j    , k]
    vzU = vz[i-1, j    , k]

    dp_dz =  (pF - pB) / ΔzC
    dσzz_dz = 2.0*(η_nF*(vzR - vzC)/ΔzR - η_nB*(vzC - vzL)/ΔzL)/ΔzC
    dσzx_dx = (
            η_xzR*((vxRF - vxLF)/ΔzC + (vzF - vzC)/ΔxF)
        - η_xzL*((vxRB - vxLB)/ΔzC + (vzC - vzB)/ΔxB)
        )/ΔxC
    dσzy_dy = (
            η_yzD*((vzD - vzC)/ΔyD + (vyLR - vyLL)/ΔzC) 
            - η_yzU*((vzC - vzU)/ΔyU + (vyUR - vyUL)/ΔzC)
            )/ΔyC

    ΔR = RZ[i,j,k] + dp_dz - dσzz_dz - dσzx_dx - dσzy_dy

    Coeff_vzC = (
        - 2.0*(η_nF/ΔzR/ΔzC +  η_nB/ΔzL/ΔzC) 
        -     (η_xzR/ΔxF/ΔxC + η_xzL/ΔxB/ΔxC) 
        -     (η_yzD/ΔyD/ΔyC + η_yzU/ΔyU/ΔyC)
        )

    return ΔR, Coeff_vzC
end

# Continuity equation dvx/dx+dvy/dy+dvz/dz=RC
# is solved via pressure updates
# dpr=-etas*div(v)
# pr-Boundary conditions 
function calculate_pressure_residual(
    i::Int64,
    j::Int64,
    k::Int64,
    level_data::LevelData
)::Float64

    vx = level_data.vx.array
    vy = level_data.vy.array
    vz = level_data.vz.array
    RC = level_data.RC.array

    xstp_b = level_data.grid.arrays.basic.xstp_b.array
    ystp_b = level_data.grid.arrays.basic.ystp_b.array
    zstp_b = level_data.grid.arrays.basic.zstp_b.array

    ΔxC = xstp_b[j]
    ΔyC = ystp_b[i]
    ΔzC = zstp_b[k]

    vxR = vx[i+1, j+1, k+1]
    vxL = vx[i+1, j  , k+1]
    vyD = vy[i+1, j+1, k+1]
    vyU = vy[i  , j+1, k+1]
    vzF = vz[i+1, j+1, k+1]
    vzB = vz[i+1, j+1, k  ]

    # Solving Continuity equation by adjusting pressure
    # Computing current residual
    ΔR = RC[i,j,k] - (
        (vxR - vxL)/ΔxC
        + (vyD - vyU)/ΔyC
        + (vzF - vzB)/ΔzC
        )

    return ΔR
end

function calculate_pressure_residual(
    i::Int64,
    j::Int64,
    level_data::LevelData2d
)::Float64

    vx = level_data.vx.array
    vy = level_data.vy.array
    RC = level_data.RC.array

    xstp_b = level_data.grid.arrays.basic.xstp_b.array
    ystp_b = level_data.grid.arrays.basic.ystp_b.array

    ΔxC = xstp_b[j]
    ΔyC = ystp_b[i]

    vxR = vx[i+1, j+1]
    vxL = vx[i+1, j  ]
    vyD = vy[i+1, j+1]
    vyU = vy[i  , j+1]

    # Solving Continuity equation by adjusting pressure
    # Computing current residual
    ΔR = RC[i,j] - (
        (vxR - vxL)/ΔxC
        + (vyD - vyU)/ΔyC
        )

    return ΔR
end

end # module