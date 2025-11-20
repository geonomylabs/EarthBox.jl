module Limits

@inline function set_yield_stress_limits(
    stress_yield_minimum::Float64,
    stress_yield_maximum::Float64,
    stress_yield::Float64
)::Float64
    stress_yield = max(stress_yield, stress_yield_minimum)
    stress_yield = min(stress_yield, stress_yield_maximum)
    return stress_yield
end

@inline function apply_viscosity_limits(
    viscosity_minimum::Float64,
    viscosity_maximum::Float64,
    viscosity::Float64
)::Float64
    viscosity = max(viscosity, viscosity_minimum)
    viscosity = min(viscosity, viscosity_maximum)
    return viscosity
end

end # module Limits 