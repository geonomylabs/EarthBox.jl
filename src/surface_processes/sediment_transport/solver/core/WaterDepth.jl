module WaterDepth

""" Calculate water depth.

# Arguments
- `topo_gridy`: Topography grid y-coordinates (meters). Note that y increases with
    depth
- `sealevel_x`: X-locations of sea level (meters)

# Returns
- `water_depth_x`: Water depth (meters)
"""
function calculate_water_depth(
    topo_gridy::Vector{Float64},
    sealevel_x::Vector{Float64}
)::Vector{Float64}
    toponum = length(topo_gridy)
    water_depth_x = zeros(toponum)
    
    for i in 1:toponum
        ytopo = topo_gridy[i]
        sealevel = sealevel_x[i]
        if sealevel < ytopo
            water_depth = ytopo - sealevel
        else
            water_depth = 0.0
        end
        water_depth_x[i] = water_depth
    end
    return water_depth_x
end

end # module 