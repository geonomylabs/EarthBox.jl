module LithosFlexure

function calc_flexural_deflection_from_block_loads(
    xcoor_blocks::Vector{Float64},
    width_blocks::Vector{Float64},
    height_blocks::Vector{Float64},
    thick_elastic_km::Float64,
    youngs_modulus::Float64,
    poissons_ratio::Float64,
    grav_acceleration::Float64,
    rho_load::Float64,
    rho_disp::Float64,
    rho_mantle::Float64,
    rho_infill::Float64
)::Vector{Float64}
    flex_rigidity = calc_flexural_rigidity(
        youngs_modulus, thick_elastic_km, poissons_ratio)
    deflections = calc_deflections_infinite_beam(
        grav_acceleration,
        flex_rigidity,
        rho_load,
        rho_disp,
        rho_mantle,
        rho_infill,
        xcoor_blocks,
        width_blocks,
        height_blocks
    )
    return deflections
end

function calc_flexural_rigidity(
    youngs_modulus::Float64,
    thick_elastic_km::Float64,
    poissons_ratio::Float64
)::Float64
    flex_rigidity = (
        youngs_modulus
        * (thick_elastic_km*1000)^3/(12*(1 - poissons_ratio^2))
    )
    return flex_rigidity
end

""" Calculate flexural deflections for an infinite beam.

Returns
-------
deflections::Vector{Float64}
    Flexural deflections in meters.
"""
function calc_deflections_infinite_beam(
    grav_acceleration::Float64,
    flex_rigidity::Float64,
    rho_load::Float64,
    rho_disp::Float64,
    rho_mantle::Float64,
    rho_infill::Float64,
    xcoor_blocks::Vector{Float64},
    width_blocks::Vector{Float64},
    height_blocks::Vector{Float64}
)::Vector{Float64}
    alpha = calc_characteristic_flexural_length_scale(
        flex_rigidity, rho_mantle, rho_infill, grav_acceleration)
    alpha_inv = 1.0/alpha

    drho1 = rho_load - rho_disp
    drho2 = rho_mantle - rho_infill
    nblocks = length(xcoor_blocks)

    deflections = zeros(nblocks)
    for i in 1:nblocks
        xloc = xcoor_blocks[i]
        deflection_sum = 0.0
        for j in 1:nblocks
            height_block = height_blocks[j]
            xblock = xcoor_blocks[j]
            half_width = width_blocks[j]/2.0
            xstart = xblock - half_width
            xend = xblock + half_width

            deflection_inc = 0.0
            if xloc <= xstart
                deflection_inc = deflection_to_left_of_block(
                    xloc, xblock, height_block, half_width,
                    drho1, drho2, alpha_inv
                )
            elseif xstart <= xloc <= xend
                deflection_inc = deflection_within_block(
                    xloc, xblock, height_block, half_width,
                    drho1, drho2, alpha_inv
                )
            else
                deflection_inc = deflection_to_right_of_block(
                    xloc, xblock, height_block, half_width,
                    drho1, drho2, alpha_inv
                )
            end
            deflection_sum = deflection_sum + deflection_inc
        end
        deflections[i] = deflection_sum
    end
    return deflections
end

function calc_characteristic_flexural_length_scale(
    flex_rigidity::Float64,
    rho_mantle::Float64,
    rho_infill::Float64,
    grav_acceleration::Float64
)::Float64
    alpha = (
        (
            4.0*flex_rigidity/(rho_mantle - rho_infill)/grav_acceleration
        )^(1.0/4.0)
    )
    return alpha
end

function deflection_to_left_of_block(
    xloc::Float64,
    xblock::Float64,
    height_block::Float64,
    half_width::Float64,
    drho1::Float64,
    drho2::Float64,
    alpha_inv::Float64
)::Float64
    deflection_inc = (
        height_block/2.0*drho1/drho2
        * (
            exp(-alpha_inv*(-xloc + xblock - half_width))
            * cos(alpha_inv*(-xloc + xblock - half_width))
            - exp(-alpha_inv*(-xloc + xblock + half_width))
            * cos(alpha_inv*(-xloc + xblock + half_width))
        )
    )
    return deflection_inc
end

function deflection_within_block(
    xloc::Float64,
    xblock::Float64,
    height_block::Float64,
    half_width::Float64,
    drho1::Float64,
    drho2::Float64,
    alpha_inv::Float64
)::Float64
    deflection_inc = (
        height_block/2.0*drho1/drho2
        * (
            2.0
            - exp(-alpha_inv*(xloc - xblock + half_width))
            * cos(alpha_inv*(xloc - xblock + half_width))
            - exp(-alpha_inv*(-xloc + xblock + half_width))
            * cos(alpha_inv*(-xloc + xblock + half_width))
        )
    )
    return deflection_inc
end

function deflection_to_right_of_block(
    xloc::Float64,
    xblock::Float64,
    height_block::Float64,
    half_width::Float64,
    drho1::Float64,
    drho2::Float64,
    alpha_inv::Float64
)::Float64
    deflection_inc = (
        -height_block/2.0*drho1/drho2
        * (
            exp(-alpha_inv*(xloc - xblock + half_width))
            * cos(alpha_inv*(xloc - xblock + half_width))
            - exp(-alpha_inv*(xloc - xblock - half_width))
            * cos(alpha_inv*(xloc - xblock - half_width))
        )
    )
    return deflection_inc
end

end # module 