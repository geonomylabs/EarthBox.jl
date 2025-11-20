module SmoothSurface

""" Smooth topography.

Inputs
------
ytopo: Array((xnum), dtype=np.float64)
    - Y-coordinates of surface.

nsmooth: int
    - Number of points to use for the running average.

Returns
-------
y_surface_smooth: Array((xnum), dtype=np.float64)
    - Smoothed y-coordinates of surface.

"""
function smooth_surface(
    ytopo::Vector{Float64};
    nsmooth::Int = 2
)::Vector{Float64}
    xnum = length(ytopo)
    y_surface_smooth = Vector{Float64}(undef, xnum)
    for i in 1:xnum
        ytopo_sum = 0.0
        icount = 0
        for j in 1:nsmooth*2 + 1
            ii = i - nsmooth + j
            if 1 <= ii <= xnum
                ytopo_sum = ytopo_sum + ytopo[ii]
                icount += 1
            end
        end
        y_surface_smooth[i] = ytopo_sum/Float64(icount)
    end
    return y_surface_smooth
end

end