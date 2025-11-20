module LayerAverage

import EarthBox: MathTools

function calculate_average(
    xtop::Float64,
    xbottom::Float64,
    xcors::Vector{Float64},
    scalars::Vector{Float64},
    dx_interp::Float64
)::Float64
    width_interp = xbottom - xtop
    nx_interp = floor(Int, width_interp/dx_interp) + 1
    if nx_interp < 0
        nx_interp = 0
        print_negative_index_warning(xtop, xbottom)
    end
    
    x_interp = zeros(nx_interp)
    scalars_interp = zeros(nx_interp)
    
    for k in 1:nx_interp
        xcor = xtop + (k-1)*dx_interp
        x_interp[k] = xcor
    end
    
    MathTools.linear_interp_vals!(xcors, scalars, x_interp, scalars_interp)
    scalar_avg = calculate_average_for_profile(x_interp, scalars_interp)
    return scalar_avg
end

function print_negative_index_warning(xtop::Int, xbottom::Int)::Nothing
    println(
        "!!! Warning !!! in calculate_average: nx_interp < 0. ",
        "xtop: ", xtop, " xbottom: ", xbottom,
        " where nx_interp is xbottom - xtop"
    )
    return nothing
end

function calculate_average_for_profile(
    xcors::Vector{Float64},
    ycors::Vector{Float64}
)::Float64
    nx = length(xcors)
    dxx = xcors[nx] - xcors[1]
    sumit = 0.0
    
    if dxx > 0.0
        for i in 2:nx
            x0 = xcors[i-1]
            x1 = xcors[i]
            dx = x1 - x0
            y0 = ycors[i-1]
            y1 = ycors[i]
            
            if y1 > y0
                area = y0*dx + (y1 - y0)*dx/2.0
            elseif y1 < y0
                area = y1*dx + (y0 - y1)*dx/2.0
            else
                area = dx*y0
            end
            sumit += area
        end
        val = sumit/dxx
    else
        val = 0.0
    end
    return val
end

end # module 