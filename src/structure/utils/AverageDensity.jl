module AverageDensity

import ..LayerAverage: calculate_average

""" Calculate average density of all layers.

Average density is calculated for vertical columns in the sticky air layer,
crustal layer, mantle lithosphere layer and asthenosphere layer for each
x-grid coordinate.
"""
function calculate_average_density_of_layers(
    xnum::Int,
    ynum::Int,
    gridy_b::Vector{Float64},
    ytopo::Vector{Float64},
    ymoho::Vector{Float64},
    ylith_base::Vector{Float64},
    rho_kg_m3::Matrix{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}, Vector{Float64}, Vector{Float64}}
    # Average density arrays
    rho_sticky_x = zeros(xnum)   # Sticky Air
    rho_crust_x = zeros(xnum)   # Crust
    rho_lithosphere_x = zeros(xnum)   # Mantle Lithosphere
    rho_asthenosphere_x = zeros(xnum)   # Asthenosphere
    
    # Vertical grid coordinates
    ycors = zeros(ynum)
    # Density at vertical grid coordinates
    density_y = zeros(ynum)
    # Y-coordinate limit used for average density of asthenosphere layer
    y_bottom = gridy_b[ynum]
    # Interpolation step size (meters)
    dy_interp = 25.0
    
    # Y-coordinate array
    for i in 1:ynum
        ycors[i] = gridy_b[i]
    end

    for j in 1:xnum
        # Get interface depths
        y_top = 0.0
        y_topo = ytopo[j]
        y_moho = ymoho[j]
        y_lith_base = ylith_base[j]
        
        # Check on lithosphere base
        if y_lith_base < y_moho
            y_lith_base = y_moho + 1000.0
        end
        
        # Define scalar values at y grid coordinates
        for i in 1:ynum
            density_y[i] = rho_kg_m3[i, j]
        end

        rho_sticky_x[j] = calculate_average(
            y_top, y_topo, ycors, density_y, dy_interp)

        rho_crust_x[j] = calculate_average(
            y_topo, y_moho, ycors, density_y, dy_interp)

        rho_lithosphere_x[j] = calculate_average(
            y_moho, y_lith_base, ycors, density_y, dy_interp)

        rho_asthenosphere_x[j] = calculate_average(
            y_lith_base, y_bottom, ycors, density_y, dy_interp)
    end

    return rho_sticky_x, rho_crust_x, rho_lithosphere_x, rho_asthenosphere_x
end

end # module 