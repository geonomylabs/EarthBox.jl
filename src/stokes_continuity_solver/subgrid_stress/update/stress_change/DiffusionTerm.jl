module DiffusionTerm

""" Calculate subgrid stress diffusion term.

# Arguments
- `viscosity`: Material viscosity
- `shear_modulus`: Material shear modulus
- `timestep`: Simulation timestep
- `subgrid_diff_coef_stress`: Subgrid diffusion coefficient for stress

# Returns
- `stress_diffusion`: Subgrid stress diffusion term
"""
@inline function calculate_subgrid_stress_diffusion_term(
    viscosity::Float64,
    shear_modulus::Float64,
    timestep::Float64,
    subgrid_diff_coef_stress::Float64
)::Float64
    # Compute local stress relaxation timescale (Maxwell time) for the marker
    maxwell_time = viscosity / shear_modulus
    # Compute degree of subgrid stress relaxation
    sdif = -subgrid_diff_coef_stress * timestep / maxwell_time
    # Limit the diffusion term to prevent numerical instabilities
    sdif = max(sdif, -30.0)
    # Calculate final diffusion term
    stress_diffusion = 1.0 - exp(sdif)
    return stress_diffusion
end

end # module 