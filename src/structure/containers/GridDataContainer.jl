module GridDataContainer

"""
    GridData

Struct used to store grid data.

# Fields
- `xnum::Union{Int, Nothing}`: Number of grid points in x-direction
- `gridx_b::Union{Vector{Float64}, Nothing}`: X-coordinates of basic grid points in meters
- `xstp_avg::Union{Float64, Nothing}`: Average basic grid spacing in x-direction in meters
- `ysize::Union{Float64, Nothing}`: Size of the model in y-direction in meters
"""
mutable struct GridData
    xnum::Union{Int, Nothing}
    gridx_b::Union{Vector{Float64}, Nothing}
    xstp_avg::Union{Float64, Nothing}
    ysize::Union{Float64, Nothing}
end

# Constructor with default values
GridData(; xnum=nothing, gridx_b=nothing, xstp_avg=nothing, ysize=nothing) = 
    GridData(xnum, gridx_b, xstp_avg, ysize)

end