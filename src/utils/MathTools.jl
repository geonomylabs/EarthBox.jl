module MathTools

import Random
import LinearAlgebra

""" Function that returns 1.0 with probability p, else 0.0.
"""
function zero_or_one(p::Float64)::Float64
    rn = Random.rand()
    return rn < p ? 1.0 : 0.0
end

""" Compute a vector vxy from vx and vy values.
"""
function compute_vxy_vec(
    vx::Matrix{Float64}, 
    vy::Matrix{Float64}, 
    xnum::Int, 
    ynum::Int, 
    vxy::Vector{Float64}
)::Nothing
    icount = 1
    # vx
    for j in 1:ynum+1
        for i in 1:xnum
            vxy[icount] = vx[j,i]
            icount += 1
        end
    end
    # vy
    for j in 1:ynum
        for i in 1:xnum+1
            vxy[icount] = vy[j,i]
            icount += 1
        end
    end
    return nothing
end

""" Calculate harmonic average of 4 values.
"""
function harmonic4(
    v1::Float64,
    v2::Float64,
    v3::Float64,
    v4::Float64,
    fac1::Float64, 
    fac2::Float64, 
    fac3::Float64, 
    fac4::Float64
)::Tuple{Float64, Int}
    vinv = 1.0/v1 + 1.0/v2 + 1.0/v3 + 1.0/v4 + fac1 + fac2 + fac3 + fac4
    icheck = 0
    if vinv != 0.0 && 4.0 + fac1 + fac2 + fac3 + fac4 != 0.0
        avg_inv = vinv / (4.0 + fac1 + fac2 + fac3 + fac4)
    else
        avg_inv = 1.0
        icheck = 1
    end
    vharm = 1.0/avg_inv
    return (vharm, icheck)
end

function get_L2_norms(
    v::Union{Vector{Float64}, Matrix{Float64}}, 
    v_old::Union{Vector{Float64}, Matrix{Float64}}
)::Tuple{Float64, Float64}
    vdiff = abs.(v_old .- v)
    vdiff_l2 = LinearAlgebra.norm(vdiff)
    v_l2 = LinearAlgebra.norm(v)
    rdiffv = v_l2 != 0.0 ? vdiff_l2 / v_l2 : 1e38
    return (vdiff_l2, rdiffv)
end

function get_inf_norms(
    v::Union{Vector{Float64}, Matrix{Float64}}, 
    v_old::Union{Vector{Float64}, Matrix{Float64}}
)::Tuple{Float64, Float64}
    vdiff = abs.(v_old .- v)
    vdiff_inf = LinearAlgebra.norm(vdiff, Inf)
    v_inf = LinearAlgebra.norm(v, Inf)
    rdiffv = v_inf != 0.0 ? vdiff_inf / v_inf : 1e38
    return (vdiff_inf, rdiffv)
end

""" Interpolate y-value at x-location given input x and y-locations.
"""
function linear_interp_at_x_location(
    x_location::Float64, 
    x_locations::Vector{Float64},
    y_locations::Vector{Float64}
)::Float64
    x_interp = [x_location]
    y_interp = zeros(1)
    linear_interp_vals!(x_locations, y_locations, x_interp, y_interp)
    return y_interp[1]
end

""" Interpolate y-values at x-locations given input x locations.
"""
function linear_interp_vals!(
    x_coors::Vector{Float64}, 
    y_values::Vector{Float64},
    x_coors_interp::Vector{Float64}, 
    y_values_interp::Vector{Float64}
)::Nothing
    nval = length(y_values)
    ninterp = length(y_values_interp)
    jmax = length(x_coors)
    ymin = x_coors[1]
    ymax = x_coors[jmax]
    for i in 1:ninterp
        yi = x_coors_interp[i]
        vi = 0.0
        for j in 2:nval
            yp_1 = x_coors[j-1]
            yp_2 = x_coors[j]
            vp_1 = y_values[j-1]
            vp_2 = y_values[j]
            if yi == yp_1
                vi = vp_1
            elseif yi == yp_2
                vi = vp_2
            elseif yi <= ymin
                vi = y_values[1]
            elseif yi >= ymax
                vi = y_values[jmax]
            elseif yp_1 < yi < yp_2
                vi = vp_1 + (vp_2 - vp_1)/(yp_2 - yp_1)*(yi - yp_1)
            end
        end
        y_values_interp[i] = vi
    end
    return nothing
end

function linear_interp_vals_opt!(
    x_coors::Vector{Float64}, 
    y_values::Vector{Float64},
    x_coors_interp::Vector{Float64}, 
    y_values_interp::Vector{Float64}
)::Nothing
    nnodes = length(y_values)
    ninterp = length(y_values_interp)
    xindex_max = length(x_coors)
    xmin = x_coors[1]
    xmax = x_coors[xindex_max]
    
    for i in 1:ninterp
        x_interpolation = x_coors_interp[i]
        value_interpolated = linear_interp(
            x_interpolation, x_coors, y_values, nnodes, xmin, xmax, xindex_max)
        y_values_interp[i] = value_interpolated
    end
    return nothing
end

""" Interpolate y-values at x-location.
"""
function linear_interp_single(
    x_coors::Vector{Float64}, 
    y_values::Vector{Float64},
    x_interpolation::Float64
)::Float64
    nnodes = length(y_values)
    xindex_max = length(x_coors)
    xmin = x_coors[1]
    xmax = x_coors[xindex_max]
    return linear_interp(x_interpolation, x_coors, y_values, nnodes, xmin, xmax, xindex_max)
end

function linear_interp(
    x_interpolation::Float64, 
    x_coors::Vector{Float64},
    y_values::Vector{Float64}, 
    nnodes::Int, xmin::Float64,
    xmax::Float64, xindex_max::Int
)::Float64
    if x_interpolation <= xmin
        return y_values[1]
    elseif x_interpolation >= xmax
        return y_values[xindex_max]
    else
        for j in 2:nnodes
            xp_1 = x_coors[j-1]
            xp_2 = x_coors[j]
            yp_1 = y_values[j-1]
            yp_2 = y_values[j]
            if x_interpolation == xp_1
                return yp_1
            elseif x_interpolation == xp_2
                return yp_2
            elseif xp_1 < x_interpolation < xp_2
                return yp_1 + (yp_2 - yp_1)/(xp_2 - xp_1)*(x_interpolation - xp_1)
            end
        end
    end
    return 0.0
end

@inline function linear_interp_bisection(
    gridx::Vector{Float64},
    gridy::Vector{Float64},
    x_location::Float64
)::Float64
    idx = searchsortedfirst(gridx, x_location) - 1
    if idx < 1
        idx = 1
    elseif idx >= length(gridx)
        idx = length(gridx) - 1
    end
    x0 = gridx[idx]
    x1 = gridx[idx + 1]
    y0 = gridy[idx]
    y1 = gridy[idx + 1]
    return x1 != x0 ? y0 + (y1 - y0) * (x_location - x0) / (x1 - x0) : y0
end

""" Generate random numbers from a normal distribution.
"""
function generate_normal_random_numbers(
    nvalues::Int; 
    mean::Float64=0.0,
    standard_deviation::Float64=1.0
)::Vector{Float64}
    return [generate_normal_random_number(mean=mean, standard_deviation=standard_deviation) for _ in 1:nvalues]
end

""" Generate a single random number from a normal distribution.
"""
function generate_normal_random_number(;
    mean::Float64=0.0,
    standard_deviation::Float64=1.0
)::Float64
    u1 = Random.rand()
    u2 = Random.rand()
    return standard_deviation * box_muller(u1, u2) + mean
end

""" Generate a random number using the Box-Muller method.
"""
function box_muller(u1::Float64, u2::Float64)::Float64
    r = sqrt(-2.0 * log(u1))
    theta = 2.0 * π * u2
    return r * sin(theta)
end

""" Smooth surface by averaging over nsmooth points.
"""
function smooth_surface(ytopo::Vector{Float64}; nsmooth::Int=2)::Vector{Float64}
    xnum = length(ytopo)
    y_surface_smooth = zeros(xnum)
    for i in 1:xnum
        ytopo_sum = 0.0
        icount = 0
        for j in -nsmooth:nsmooth
            ii = i + j
            if 1 <= ii <= xnum
                ytopo_sum += ytopo[ii]
                icount += 1
            end
        end
        y_surface_smooth[i] = ytopo_sum/icount
    end
    return y_surface_smooth
end

end # module 