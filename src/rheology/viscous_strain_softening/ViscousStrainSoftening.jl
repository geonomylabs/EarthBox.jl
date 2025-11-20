module ViscousStrainSoftening

"""
    update_pre_exponential_for_viscous_strain_softening(
        pre_exponential_dislocation::Float64,
        strain::Float64,
        vsoftfac::Float64,
        strain_initial::Float64,
        strain_final::Float64
    )::Float64

Update pre-exponential term for power-law due to strain softening.
"""
@inline function update_pre_exponential_for_viscous_strain_softening(
    pre_exponential_dislocation::Float64,
    strain::Float64,
    vsoftfac::Float64,
    strain_initial::Float64,
    strain_final::Float64
)::Float64
    pre_exponential_dislocation_min = pre_exponential_dislocation
    pre_exponential_dislocation_max = pre_exponential_dislocation * vsoftfac
    
    if strain_initial < strain < strain_final
        pre_exponential_dislocation = (
            pre_exponential_dislocation_min
            + (
                pre_exponential_dislocation_max
                - pre_exponential_dislocation_min
            ) / (strain_final - strain_initial) * (strain - strain_initial)
        )
    end
    
    if strain >= strain_final
        pre_exponential_dislocation = pre_exponential_dislocation_max
    end
    
    return pre_exponential_dislocation
end

end # module ViscousStrainSoftening 