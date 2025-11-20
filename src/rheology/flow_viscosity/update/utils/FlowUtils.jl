module FlowUtils

"""
    apply_minimum_temperature_limit(temperature::Float64, temperature_minimum::Float64)::Float64

Limit marker temperature to avoid numerical issues with flow laws.
"""
function apply_minimum_temperature_limit(
    temperature::Float64,
    temperature_minimum::Float64
)::Float64
    return max(temperature, temperature_minimum)
end


end # module

