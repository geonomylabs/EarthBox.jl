module WaterDepth

""" Calculate water depth.

# Arguments
- `topo_gridy`: Topography grid y-coordinates (meters). Note that y increases with
    depth
- `sealevel_x`: X-locations of sea level (meters)

# Returns
- `water_depth_x`: Water depth (meters)
"""
function calculate_water_depth!(
    water_depth_x::Vector{Float64},
    topo_gridy::Vector{Float64},
    sealevel_x::Vector{Float64}
)::Nothing
    toponum = length(topo_gridy)
    for i in 1:toponum
        ytopo = topo_gridy[i]
        sealevel = sealevel_x[i]
        if sealevel < ytopo
            water_depth_x[i] = ytopo - sealevel
        else
            water_depth_x[i] = 0.0
        end
    end
    return nothing
end

end # module 