module Diffusivity

import EarthBox.ConversionFuncs: mm_per_yr_to_meters_per_seconds
import EarthBox.DataStructures: SedimentTransportParameters

function make_topo_diffusivity_grid(
    topo_gridx::Vector{Float64},
    water_depth_x::Vector{Float64},
    downstream_distances_x::Vector{Float64},
    sediment_transport_parameters::SedimentTransportParameters,
    use_constant_diffusivity::Bool
)::Vector{Float64}
    if use_constant_diffusivity
        xnum = length(topo_gridx)
        topo_diff_coef = sediment_transport_parameters.subaerial_slope_diffusivity
        return fill(topo_diff_coef, xnum)
    else
        return calculate_sediment_transport_diffusivity(
            topo_gridx, water_depth_x,
            downstream_distances_x,
            sediment_transport_parameters.subaerial_slope_diffusivity,
            sediment_transport_parameters.precipitation_rate,
            sediment_transport_parameters.subaerial_transport_coefficient,
            sediment_transport_parameters.submarine_slope_diffusivity,
            sediment_transport_parameters.submarine_diffusion_decay_depth
        )
    end
end

""" Calculate slope diffusion coefficient.

Typical values for topography diffusion coefficient range from 0 to 10^5 m^2/year.

# Arguments
- `transport_length_scale`: Transport length scale (meters)
- `maximum_elevation`: Maximum elevation (meters)
- `erosion_velocity_mm_yr`: Erosion velocity (mm/year)

# Returns
- `topo_diff_coef`: Topography diffusion coefficient (m^2/s)
"""
function calculate_slope_diffusion_coefficient(
    transport_length_scale::Float64,
    maximum_elevation::Float64,
    erosion_velocity_mm_yr::Float64
)::Float64
    erosion_velocity_m_s = mm_per_yr_to_meters_per_seconds(erosion_velocity_mm_yr)
    return erosion_velocity_m_s * transport_length_scale^2.0 / maximum_elevation
end

""" Calculate sediment transport diffusivity.

# Arguments
- `topo_gridx`: Topography grid x-coordinates (meters)
- `water_depth_x`: Water depth (meters)
- `downstream_distances`: Downstream distances for each node (meters). This 
    parameter represents the width of the drainage basin
- `subaerial_slope_diffusivity`: Subaerial slope diffusivity (m^2/s)
- `precipitation_rate`: Precipitation rate (m/s)
- `subaerial_transport_coefficient`: Subaerial transport coefficient
- `submarine_slope_diffusivity`: Submarine slope diffusivity (m^2/s)
- `submarine_diffusion_decay_depth`: Submarine diffusion decay depth (meters)

# Returns
- `topo_grid_diffusivity`: Topography diffusion coefficient (m^2/s)
"""
function calculate_sediment_transport_diffusivity(
    topo_gridx::Vector{Float64},
    water_depth_x::Vector{Float64},
    downstream_distances::Vector{Float64},
    subaerial_slope_diffusivity::Float64,
    precipitation_rate::Float64,
    subaerial_transport_coefficient::Float64,
    submarine_slope_diffusivity::Float64,
    submarine_diffusion_decay_depth::Float64
)::Vector{Float64}
    toponum = length(topo_gridx)
    topo_grid_diffusivity = zeros(toponum)
    
    for i in 2:toponum
        dd_dist = downstream_distances[i]
        water_depth = water_depth_x[i]
        if water_depth > 0.0
            topo_grid_diffusivity[i] = submarine_slope_diffusivity * 
                exp(-water_depth / submarine_diffusion_decay_depth)
        else
            topo_grid_diffusivity[i] = subaerial_slope_diffusivity + 
                precipitation_rate * subaerial_transport_coefficient * dd_dist
        end
    end
    return topo_grid_diffusivity
end

""" Print topography diffusion coefficient in mm^2/yr.
"""
function print_diffusion_coefficient(topo_diff_coef::Float64)::Nothing
    println(
        "Topography diffusion coefficient m^2/yr: ",
        topo_diff_coef * (365 * 24 * 3600)
    )
    return nothing
end

end # module 