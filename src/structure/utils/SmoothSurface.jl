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
    smooth_surface!(y_surface_smooth, ytopo, nsmooth=nsmooth)
    return y_surface_smooth
end

""" Smooth topography in place, writing to a pre-allocated output buffer.

The output buffer must be a different array from the input (the algorithm
reads `ytopo[i ± nsmooth]` while writing `out[i]`; aliasing would corrupt
results). Both arrays must have the same length.

# Arguments
- `out::Vector{Float64}`: Pre-allocated output buffer, same length as `ytopo`.
- `ytopo::Vector{Float64}`: Input surface y-coordinates.
- `nsmooth::Int`: Number of points on each side of `i` to use in the running
  average.
"""
function smooth_surface!(
    out::Vector{Float64},
    ytopo::Vector{Float64};
    nsmooth::Int = 2
)::Nothing
    @assert out !== ytopo "smooth_surface! requires distinct out and ytopo"
    xnum = length(ytopo)
    @assert length(out) == xnum "smooth_surface! requires length(out) == length(ytopo)"
    @inbounds for i in 1:xnum
        ytopo_sum = 0.0
        icount = 0
        for j in 1:nsmooth*2 + 1
            ii = i - nsmooth + j
            if 1 <= ii <= xnum
                ytopo_sum = ytopo_sum + ytopo[ii]
                icount += 1
            end
        end
        out[i] = ytopo_sum/Float64(icount)
    end
    return nothing
end

end