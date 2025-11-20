module MeltPropertiesOpt

import ..MeltFraction

""" Calculate melting-related rock properties for marker.

This function is used to include the effects of melting on density (including 
partially molten mixtures and depletion of the matrix), and the effects of 
latent heating with temperature and pressure changes over time.

# Arguments
- `pressure_pascals::Float64`: Pressure in Pascals
- `temperature_kelvins::Float64`: Temperature in Kelvins
- `density::Float64`: Density of marker (kg/m^3)
- `rhocp::Float64`: Product of marker density (kg/m^3) and heat capacity (J/K/kg)
- `adiabatic_term::Float64`: Marker adiabatic heating term (expansivity x temp)
- `iuse_melt_thermal_props::Int64`: Flag to use melt thermal properties
- `density_melt::Float64`: Density of melt (kg/m^3)
- `_melt_fraction::Float64`: Melt fraction
- `extracted_melt_fraction::Float64`: Extracted melt fraction
- `extractable_melt_fraction::Float64`: Extractable melt fraction
- `itype_solidus::Int`: Solidus model index
- `itype_liquidus::Int`: Liquidus model index
- `latent_heat::Float64`: Latent heat of melting (J/kg)
- `iuse_depletion_density::Int64`: Flag to use depletion density

# Returns
- `density::Float64`: Updated density of marker (kg/m^3)
- `rhocp::Float64`: Updated product of marker density and heat capacity
- `adiabatic_term::Float64`: Updated marker adiabatic heating term
- `delta_heat_capacity::Float64`: Change in effective heat capacity
- `delta_expansivity::Float64`: Change in effective expansivity
"""
@inline function marker_melt_props(
    pressure_pascals::Float64,
    temperature_kelvins::Float64,
    density::Float64,
    rhocp::Float64,
    adiabatic_term::Float64,
    iuse_melt_thermal_props::Int64,
    density_melt::Float64,
    _melt_fraction::Float64,
    extracted_melt_fraction::Float64,
    extractable_melt_fraction::Float64,
    itype_solidus::Int,
    itype_liquidus::Int,
    latent_heat::Float64,
    iuse_depletion_density::Int64
)::Tuple{Float64, Float64, Float64, Float64, Float64}
    if iuse_depletion_density == 1
        density += calculated_delta_density_for_depletion(extracted_melt_fraction)
    end

    delta_heat_capacity = 0.0
    delta_expansivity = 0.0
    if iuse_melt_thermal_props == 1
        if melt_and_extractable(_melt_fraction, extractable_melt_fraction)

            delta_meltfrac = _melt_fraction - extracted_melt_fraction

            density, rhocp = update_bulk_properties_based_on_melt_solid_mixing(
                density,
                rhocp,
                density_melt,
                pressure_pascals,
                delta_meltfrac
            )

            if latent_heating_potential(_melt_fraction)
                adiabatic_term, delta_expansivity = 
                    update_adiabatic_heating_term_to_included_latent_heating(
                        adiabatic_term,
                        density,
                        pressure_pascals,
                        temperature_kelvins,
                        itype_solidus,
                        itype_liquidus,
                        latent_heat
                    )

                rhocp, delta_heat_capacity = update_rhocp_to_include_latent_heating(
                    pressure_pascals,
                    temperature_kelvins,
                    rhocp,
                    density,
                    itype_solidus,
                    itype_liquidus,
                    latent_heat
                )
            end
        end
    end

    return density, rhocp, adiabatic_term, delta_heat_capacity, delta_expansivity
end

""" Check if melt fraction and extractable melt > 0.

Only use melting model if melt fraction is greater than 0 and incremental
melt extraction >= 0. Otherwise the material is refractory.

# Arguments
- `_melt_fraction::Float64`: Melt fraction
- `extractable_melt_fraction::Float64`: Extractable melt fraction

# Returns
- `check::Bool`: True if melt fraction > 0 and extractable melt fraction >= 0
"""
@inline function melt_and_extractable(
    _melt_fraction::Float64,
    extractable_melt_fraction::Float64
)::Bool
    return _melt_fraction > 0 && extractable_melt_fraction >= 0
end

""" Check for the potential for latent heating.

# Arguments
- `_melt_fraction::Float64`: Melt fraction

# Returns
- `check::Bool`: True if melt fraction < 1.0
"""
@inline function latent_heating_potential(_melt_fraction::Float64)::Bool
    return _melt_fraction < 1.0
end

""" Update density and rho*cp for melting.

# Arguments
- `density::Float64`: Density of marker (kg/m^3)
- `rhocp::Float64`: Product of marker density and heat capacity
- `density_melt::Float64`: Density of melt (kg/m^3)
- `pressure_pascals::Float64`: Pressure in Pascals
- `delta_meltfrac::Float64`: Change in melt fraction

# Returns
- `density::Float64`: Updated density of marker (kg/m^3)
- `rhocp::Float64`: Updated product of marker density and heat capacity
"""
function update_bulk_properties_based_on_melt_solid_mixing(
    density::Float64,
    rhocp::Float64,
    density_melt::Float64,
    pressure_pascals::Float64,
    delta_meltfrac::Float64
)::Tuple{Float64, Float64}
    density_melt_eos = calculate_melt_density_birch_murnagham(
        density_melt, pressure_pascals)
    density = calculate_property_for_rock_melt_mixture(
        delta_meltfrac, density_melt_eos, density)
    rhocp = calculate_property_for_rock_melt_mixture(
        delta_meltfrac, density_melt_eos, rhocp)
    return density, rhocp
end

""" Linearized Birch-Murnagham EOS for Anhydrous Basalt

The melt density equation is as follows:

    density_melt_eos = density_melt_1bar + pressure_pascals*f

where f is 8.97e-8 kg/m^3/Pa. The target pressure range is associated with
depths of 50-55 km. Constant was derived using Kto = 20.8, Ktp = 4.6
(Lesher and Spera, 2015, Table S5.4; See also Stopler et al. 1981).

Equation 1 also provides a reasonable approximation for any silicate melt.

# Arguments
- `density_melt::Float64`: Density of melt at 1 bar (kg/m^3)
- `pressure_pascals::Float64`: Pressure in Pascals

# Returns
- `density_melt_eos::Float64`: Updated density of melt (kg/m^3)

# References
- Lesher, C. E., Spera, F. J. (2015), Thermodynamic and transport properties of
  silicate melts and magma. in The Encyclopedia of Volcanoes, 2nd Edition.
- Stopler, E. et al. (1981). Melt Segregation from partially molten source
  regions: the importance of melt density and source region size. JGR, 86,
  6261-6271.
"""
@inline function calculate_melt_density_birch_murnagham(
    density_melt::Float64,
    pressure_pascals::Float64
)::Float64
    delta_pressure_eos = 1.5e9
    # density increase (kg/m3) of silicate melt at 1.5 Ga
    delta_density_eos = 130.0
    density_melt_eos = density_melt + pressure_pascals * delta_density_eos / 
        delta_pressure_eos
    return density_melt_eos
end

""" Calculate new property based on melt fraction.

# Arguments
- `delta_meltfrac::Float64`: Change in melt fraction
- `property_melt::Float64`: Property of melt
- `property_solid::Float64`: Property of solid

# Returns
- `property_mixture::Float64`: Updated property
"""
@inline function calculate_property_for_rock_melt_mixture(
    delta_meltfrac::Float64,
    property_melt::Float64,
    property_solid::Float64
)::Float64
    property_mixture = property_solid * (1.0 - delta_meltfrac) + 
        property_melt * delta_meltfrac
    return property_mixture
end

""" Update adiabatic heat term based on latent heat.

A simple 2D heat heat equation with conduction, adiabatic heating/cooling
and latent heat is given by:

    rho*Cp*dT/dt = k*laplacian(T) + Ha + Hm                            eq. 1

where rho is density, Cp is heat capacity, k is thermal conductivity, T is
temperature, Hm is the latent heat source/sink term associated with
crystallization and melting and Ha is the adiabatic heat source/sink term. 

The adiabatic term Ha is given by:

    Ha = alpha*T*(dP/dt)                                               eq. 2

where alpha is thermal expansivity, T is temperature, and dP/dt is the rate
of change of pressure with time. Note if dP/dt is positive, then the 
adiabatic term is positive and the material undergoes adiabatic heating.
If dP/dt is negative, then the adiabatic term is negative and the material
undergoes adiabatic cooling. Note also the units of equation 2 are W/m^3.

The latent heating/cooling term Hm is given by:

    Hm = -rho*L*(dM/dt)                                                eq. 3

where L is latent heat, dM/dt is the rate of change of melt fraction with
time. If dM/dt is negative (i.e. melt fraction is decreasing), then latent
heat is released due to crystallization and the heat source term is positive.
If dM/dt is positive (i.e. melt fraction is increasing), then latent heat is
absorbed due to melting and the heat source term is negative. Note that the 
units of equation 2 are W/m^3.

Equation 3 to can be reformulated as:

    Hm = -rho*L*(dM/dP)*(dP/dt)                                        eq. 3

where dM/dP is the derivative of melt fraction with respect to pressure at
constant temperature. Note that the units of equation 3 are still W/m^3.

Equation 2 and 3 can be used to rewrite Equation 1 as:

    rho*Cp*dT/dt = k*laplacian(T) + alpha*T*(dP/dt) - rho*L*(dM/dP)*(dP/dt) eq. 4

which in turn can be rewritten as:

    rho*Cp*dT/dt = k*laplacian(T) + [alpha*T - rho*L*(dM/dP)]*(dP/dt)  eq. 5

It is interesting to note that the sign in front of the latent heat term
in equations 4 and 5 is negative. This is because if the melt fraction
increases with decreasing pressure (i.e. dM/dP < 0) and pressure is
decreasing with time (i.e. dP/dt < 0), the heat is absorbed due to latent
heat of melting and the material cools (i.e. a heat sink develops).

Note also how the combined adiabatic and latent heat terms in equation 5 can
be expressed as an effective thermal expansivity alpha_eff:

    alpha_eff = alpha - rho*L*(dM/dP)/T                                eq. 6

allowing equation 5 to be rewritten as:

    rho*Cp*dT/dt = k*laplacian(T) + alpha_eff*T*(dP/dt)                eq. 7

We note that equation 6 differs from the formulation in Gerya (2019) where

    alpha_eff = alpha + rho*L*(dM/dP)/T                                eq. 8

We view equation 6 as the correct formulation since it is consistent with
the sign convention for latent heat in equation 3.

We also note that in the MatLab code of Gerya (2019) the following 
formulation is used:

    % Melting adiabatic term: alpham=-rho*(dHlat/dP)/T
    % Numerical differentiation
    dp=1000; % Pressure increment, Pa
    [xmelt hlat0]=Melt_fraction(MPR(mm1)-dp,MTK(mm1),MI(mm1));
    [xmelt hlat1]=Melt_fraction(MPR(mm1)+dp,MTK(mm1),MI(mm1));
    MHACUR=MHACUR-(hlat1-hlat0)/(2.0*dp);
    
This formulation for the adiabatic term with latent heating does not include
density as is properly treated in the MATLAB comment and, therefore,
underestimates the effect of latent heat.

"""
function update_adiabatic_heating_term_to_included_latent_heating(
    adiabatic_term::Float64,
    density::Float64,
    pressure_pascals::Float64,
    temperature_kelvins::Float64,
    itype_solidus::Int,
    itype_liquidus::Int,
    latent_heat::Float64
)::Tuple{Float64, Float64}
    delta_pressure = 1000.0  # Pa

    # Note that the change in melt fraction is calculated over 2*delta_pressure
    delta_xmelt_times_latent_heat = 
        calculate_delta_meltfrac_times_latent_heat_for_pressure(
            delta_pressure,
            pressure_pascals,
            temperature_kelvins,
            itype_solidus,
            itype_liquidus,
            latent_heat
        )

    # rho*Cp*dT/dt = k*laplacian(T) + [alpha*T - rho*L*(dM/dP)]*(dP/dt)
    latent_heating_term = -density * (delta_xmelt_times_latent_heat) / 
        (2.0 * delta_pressure)

    # Note that the adiabatic heating term will be multiplied by dP/dt
    # when constructing the right-hand side of the heat equation
    adiabatic_term += latent_heating_term

    # This is the effect expansivity change associated with latent heat
    delta_expansivity = latent_heating_term / temperature_kelvins

    return adiabatic_term, delta_expansivity
end

""" Calculate the change in melt fraction associated with pressure change.

The change in melt fraction is calculated over 2*delta_pressure.

Latent heat is multiplied by the change in melt fraction for use in
equations that describe how melting affects thermal expansivity.

# Arguments
- `delta_pressure::Float64`: Change in pressure (Pa)
- `pressure_pascals::Float64`: Pressure in Pascals
- `temperature_kelvins::Float64`: Temperature in Kelvins
- `itype_solidus::Int`: Solidus model index
- `itype_liquidus::Int`: Liquidus model index
- `latent_heat::Float64`: Latent heat of melting (J/kg)

# Returns
- `delta_xmelt_times_latent_heat::Float64`: Change in xmelt times latent heat
"""
function calculate_delta_meltfrac_times_latent_heat_for_pressure(
    delta_pressure::Float64,
    pressure_pascals::Float64,
    temperature_kelvins::Float64,
    itype_solidus::Int,
    itype_liquidus::Int,
    latent_heat::Float64
)::Float64
    _xmelt, xmelt_times_latent_heat_0 = 
        MeltFraction.calc_melt_fraction_and_melt_fraction_times_latent_heat(
            pressure_pascals - delta_pressure,
            temperature_kelvins,
            itype_solidus,
            itype_liquidus,
            latent_heat
        )

    _xmelt, xmelt_times_latent_heat_1 = 
        MeltFraction.calc_melt_fraction_and_melt_fraction_times_latent_heat(
            pressure_pascals + delta_pressure,
            temperature_kelvins,
            itype_solidus,
            itype_liquidus,
            latent_heat
        )

    delta_melt_fraction_times_latent_heat = 
        xmelt_times_latent_heat_1 - xmelt_times_latent_heat_0
    return delta_melt_fraction_times_latent_heat
end

""" Update rhocp based on latent heat.

The 2D heat conduction equation with latent heating can be written as:

    rho*Cp*dT/dt = k*laplacian(T) - rho*L*(dM/dt)                      eq. 1

where rho is density, Cp is heat capacity, k is thermal conductivity, T is
temperature, L is latent heat, and dM/dt is the rate of change of melt
fraction with time. Note that the units of the latent heat term are W/m^3.

The latent heat term can be rewritten as:

    -rho*L*(dM/dt) = -rho*L*(dM/dT)*(dT/dt)                            eq. 2

where dM/dT is the derivative of melt fraction with respect to temperature
at constant pressure. Note that the units of the latent heat term are still
W/m^3.

Equation 2 can be used to rewrite Equation 1 as:

    rho*Cp*dT/dt = k*laplacian(T) - rho*L*(dM/dT)*(dT/dt)              eq. 3

which can be rewritten as:

    rho*Cp*dT/dt + rho*L*(dM/dT)*dT/dt = k*laplacian(T)                eq. 4

    [Cp + L*(dM/dT)]*rho*dT/dt = k*laplacian(T)                        eq. 5

Thus the effect of latent heating with temperature changes over time can be
expressed as an effective heat capacity Cp_eff:

    Cp_eff = Cp + L*(dM/dT)                                            eq. 6

In terms of an rho*cp term, the effect of latent heating with temperature
changes over time can be expressed as:

rho*Cp_eff = rho*Cp + rho*L*(dM/dT)                                    eq. 7

Equation 7 is identical to the one used by the Matlab code of Gerya (2019).

Inputs
------
pressure_pascals: float
    Pressure in Pascals.

temperature_kelvins: float
    Temperature in Kelvins.

rhocp: float
    Initial product of marker density (kg/m^3) and heat capacity (J/K/kg).

density: float
    Density of marker (kg/m^3).

itype_solidus: int
    Solidus model index.

itype_liquidus: int
    Liquidus model index.

latent_heat: float
    Latent heat of melting (J/kg).

Returns
-------
rhocp: float
    Updated product of marker density (kg/m^3) and heat capacity (J/K/kg).

delta_heat_capacity: float
    Change in heat capacity associated with latent heat.

"""
function update_rhocp_to_include_latent_heating(
    pressure_pascals::Float64,
    temperature_kelvins::Float64,
    rhocp::Float64,
    density::Float64,
    itype_solidus::Int,
    itype_liquidus::Int,
    latent_heat::Float64
)::Tuple{Float64, Float64}
    delta_temperature = 1.0  # Kelvins

    # Note that the delta xmelt is computed over 2*delta_temperature
    delta_xmelt_times_latent_heat = 
        calculate_delta_meltfrac_times_latent_heat_for_temperature(
            pressure_pascals,
            temperature_kelvins,
            delta_temperature,
            itype_solidus,
            itype_liquidus,
            latent_heat
        )

    delta_heat_capacity = (delta_xmelt_times_latent_heat) / 
        (2.0 * delta_temperature)
    rhocp += density * delta_heat_capacity
    return rhocp, delta_heat_capacity
end

""" Calculate the change in xmelt associated with temperature change.

Latent heat is multiplied by the change in melt fraction for use in
equations that describe how melting affects heat capacity.

# Arguments
- `pressure_pascals::Float64`: Pressure in Pascals
- `temperature_kelvins::Float64`: Temperature in Kelvins
- `delta_temperature::Float64`: Change in temperature (K)
- `itype_solidus::Int`: Solidus model index
- `itype_liquidus::Int`: Liquidus model index
- `latent_heat::Float64`: Latent heat of melting (J/kg)

# Returns
- `delta_xmelt_times_latent_heat::Float64`: Change in xmelt times latent heat
"""
function calculate_delta_meltfrac_times_latent_heat_for_temperature(
    pressure_pascals::Float64,
    temperature_kelvins::Float64,
    delta_temperature::Float64,
    itype_solidus::Int,
    itype_liquidus::Int,
    latent_heat::Float64
)::Float64
    _xmelt, xmelt_times_latent_heat_0 = 
        MeltFraction.calc_melt_fraction_and_melt_fraction_times_latent_heat(
            pressure_pascals,
            temperature_kelvins - delta_temperature,
            itype_solidus,
            itype_liquidus,
            latent_heat
        )

    _xmelt, xmelt_times_latent_heat_1 = 
        MeltFraction.calc_melt_fraction_and_melt_fraction_times_latent_heat(
            pressure_pascals,
            temperature_kelvins + delta_temperature,
            itype_solidus,
            itype_liquidus,
            latent_heat
        )

    delta_xmelt_times_latent_heat = 
        xmelt_times_latent_heat_1 - xmelt_times_latent_heat_0
    return delta_xmelt_times_latent_heat
end

""" Calculate the change in density associated with melt extraction.

See Theunissen et al. (2022) for details.

# Arguments
- `extracted_melt_fraction::Float64`: Extracted melt fraction

# Returns
- `delta_density::Float64`: Change in density associated with melt extraction
"""
function calculated_delta_density_for_depletion(
    extracted_melt_fraction::Float64
)::Float64
    if extracted_melt_fraction < 0.075
        delta_density = -3.8 * extracted_melt_fraction * 100.0
    else
        delta_density = -(24.84 + 0.488 * extracted_melt_fraction * 100.0)
    end
    return delta_density
end

end # module 