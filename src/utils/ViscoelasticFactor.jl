module ViscoelasticFactor

function calc_viscoelastic_factor(
    i::Int,
    j::Int,
    viscosity::Array{Float64,2},
    shear_modulus::Array{Float64,2},
    timestep::Float64
)::Float64
    return calculate_viscoelastic_factor(viscosity[i,j], shear_modulus[i,j], timestep)
end

""" Calculate viscoelastic factor used in Stokes equations.

See pages 179 and 182-183 of Gerya (2010) for a detailed description of
the viscoelastic_factor. Note that viscoelastic_factor = (1 - Z) where Z is
defined in equation 13.3 on page 179 of Gerya (2010).
"""
function calculate_viscoelastic_factor(
    viscosity::Float64,
    shear_modulus::Float64,
    timestep::Float64
)::Float64
    denominator = viscosity + timestep * shear_modulus
    if denominator == 0
        print_warning(viscosity, shear_modulus, timestep)
    end
    return viscosity / denominator
end

function print_warning(
    viscosity::Float64,
    shear_modulus::Float64,
    timestep::Float64
)::Nothing
    println("!!! Warning !!! Denominator is zero in for viscoelastic factor.")
    println("viscosity = ", viscosity)
    println("shear_modulus = ", shear_modulus)
    println("timestep = ", timestep)
    return nothing
end

end # module 