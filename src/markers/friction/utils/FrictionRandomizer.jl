module FrictionRandomizer

"""
    randomize_initial_friction_coefficient(
        friction_coefficient::Float64,
        delta_fric_coef::Float64,
        random_number::Float64
    )::Float64

Randomize initial friction coefficient.
"""
@inline function randomize_initial_friction_coefficient(
    friction_coefficient::Float64,
    delta_fric_coef::Float64, 
    random_number::Float64
)::Float64
    if random_number <= 0.5
        friction_coefficient = friction_coefficient - delta_fric_coef
    end
    return friction_coefficient
end

end # module 