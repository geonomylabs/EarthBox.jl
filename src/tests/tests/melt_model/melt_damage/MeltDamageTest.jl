""" Test melt damage model

In the main code, Cohesion and friction coefficients are divided by the melt
damage factor to account for the reduction in strength due to melt damage.
A melt damage of 1 means no damage and a melt damage greater than 1 means
damage.
"""
module MeltDamageTest

using Plots
import EarthBox.MeltModel.MeltDamage: calculate_damage_factor_cos
import EarthBox.MeltModel.MeltDamage: calculate_damage_factor_probabilistic
import EarthBox.MeltModel.MeltDamage: calculate_damage_probability
import EarthBox.MeltModel.MeltDamage: linear_melt_damage_probability_model

function run_test()::Nothing
    avg_shallow_partial_melt_xcoor_mantle = 250_000.0
    melt_damage_distance = 20_000.0
    melt_damage_factor = 10.0
    maximum_damage_probability = 0.8

    xo_plot = avg_shallow_partial_melt_xcoor_mantle - 40_000.0
    xf_plot = avg_shallow_partial_melt_xcoor_mantle + 40_000.0

    xo = 0.0
    xf = 500_000.0
    dx = 500.0
    nx = floor(Int, (xf - xo) / dx) + 1

    xcoors = range(xo, xf, length=nx)
    damage_array = ones(nx)
    prob_array = zeros(nx)

    use_prob = true
    for i in 1:nx
        x_marker = xcoors[i]
        if use_prob
            damage_array[i] = calculate_damage_factor_probabilistic(
                x_marker, avg_shallow_partial_melt_xcoor_mantle,
                melt_damage_distance, melt_damage_factor,
                maximum_damage_probability
            )
            prob_array[i] = calculate_damage_probability(
                x_marker, avg_shallow_partial_melt_xcoor_mantle,
                melt_damage_distance, maximum_damage_probability
            )
        else
            damage_array[i] = calculate_damage_factor_cos(
                x_marker, avg_shallow_partial_melt_xcoor_mantle,
                melt_damage_distance, melt_damage_factor
            )
        end
    end

    p1 = plot(xcoors, damage_array, xlims=(xo_plot, xf_plot), ylims=(0.0, 15.0))
    savefig(p1, "melt_damage.png")

    p2 = plot(xcoors, prob_array, xlims=(xo_plot, xf_plot), ylims=(0.0, 1.0))
    savefig(p2, "melt_damage_prob.png")

    magmatic_crust_height_threshold = 500.0
    magmatic_crust_height_minimum = 750.0
    magmatic_crust_height_intermediate = 2_000.0
    magmatic_crust_height_maximum = 3_000.0
    intermediate_damage_probability = 0.1
    maximum_damage_probability = 0.8

    xtho = 0.0
    xthf = 15_000.0
    xth_pts = range(xtho, xthf, length=500)
    npts = length(xth_pts)
    prob_array = zeros(npts)
    
    for i in 1:npts
        magmatic_crust_height = xth_pts[i]
        prob = linear_melt_damage_probability_model(
            magmatic_crust_height,
            magmatic_crust_height_threshold,
            magmatic_crust_height_minimum,
            magmatic_crust_height_intermediate,
            magmatic_crust_height_maximum,
            maximum_damage_probability,
            intermediate_damage_probability
        )
        prob_array[i] = prob
    end
    
    p3 = plot(xth_pts, prob_array, xlims=(xtho, xthf), ylims=(0.0, 1.0))
    savefig(p3, "melt_damage_prob_linear.png")
    
    return nothing
end

end # module 