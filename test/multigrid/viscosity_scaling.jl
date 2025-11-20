module ViscosityScalingTest

function calculate_scaled_viscosity()
    η_min_cur = 1e20
    η_max_cur = 1e23
    η_min = 1e20
    η_max = 1e23
    η_cur = 1e21

    η_coeff = log(η_max_cur / η_min_cur) / log(η_max / η_min)
    println("η_coeff = $η_coeff")
    
    factor = exp.(log.(η_cur / η_min) * η_coeff)
    println("factor = $factor")
    
    η_scaled = η_min_cur * factor
    return η_scaled
end

end