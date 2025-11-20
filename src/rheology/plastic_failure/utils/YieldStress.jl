module YieldStress

using ..FluidPressure

"""
    calc_yield_stress(
        cohesion::Float64,
        sine_of_friction_angle::Float64,
        pressure::Float64,
        iuse_fluid_pressure_for_yield::Int,
        y_location::Float64,
        y_sealevel::Float64
    )::Float64

Calculate yield stress.

# Arguments
- `cohesion::Float64`: Cohesion in Pa
- `sine_of_friction_angle::Float64`: Sine of friction angle
- `pressure::Float64`: Pressure in Pa
- `iuse_fluid_pressure_for_yield::Int`: Flag to use fluid pressure (0 = no, 1 = yes)
- `y_location::Float64`: Y-coordinate of location in m
- `y_sealevel::Float64`: Y-coordinate of sealevel in m

# Returns
- `stress_yield::Float64`: Yield stress in Pa
"""
@inline function calc_yield_stress(
    cohesion::Float64,
    sine_of_friction_angle::Float64,
    pressure::Float64,
    iuse_fluid_pressure_for_yield::Int,
    y_location::Float64,
    y_sealevel::Float64
)::Float64
    friction_angle_radians = asin(sine_of_friction_angle)
    cosine_friction_angle = cos(friction_angle_radians)
    
    if iuse_fluid_pressure_for_yield == 0
        stress_yield = (
            cohesion * cosine_friction_angle
            + sine_of_friction_angle * pressure
        )
    else
        pressure_fluid = FluidPressure.calculate_fluid_pressure(y_location, y_sealevel)
        if pressure < pressure_fluid
            use_method1 = false
            if use_method1
                stress_yield = stress_yield_tension_method1(
                    cohesion, cosine_friction_angle,
                    pressure, pressure_fluid
                )
            else
                stress_yield = stress_yield_tension_method2(
                    cohesion, cosine_friction_angle
                )
            end
        else
            stress_yield = (
                cohesion * cosine_friction_angle
                + sine_of_friction_angle * (pressure - pressure_fluid)
            )
        end
    end
    stress_yield = max(stress_yield, 0.0)
    return stress_yield
end

"""
    stress_yield_tension_method1(
        cohesion::Float64,
        cosine_friction_angle::Float64,
        pressure::Float64,
        pressure_fluid::Float64
    )::Float64

Calculate yield stress under tension using method 1.

See Gerya 2013.
"""
@inline function stress_yield_tension_method1(
    cohesion::Float64,
    cosine_friction_angle::Float64,
    pressure::Float64,
    pressure_fluid::Float64
)::Float64
    stress_yield = (
        cohesion * cosine_friction_angle
        + 1.0 * (pressure - pressure_fluid)
    )
    return stress_yield
end

"""
    stress_yield_tension_method2(
        cohesion::Float64,
        cosine_friction_angle::Float64
    )::Float64

Calculate yield stress under tension using method 2.
"""
@inline function stress_yield_tension_method2(
    cohesion::Float64,
    cosine_friction_angle::Float64
)::Float64
    return cohesion * cosine_friction_angle
end

end # module 