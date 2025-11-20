module ArrayCollection

include("group/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractArrayCollection
import .RhsHeatGroup: RhsHeat
import .RadiogenicProductionGroup: RadiogenicProduction
import .AdiabaticProductionGroup: AdiabaticProduction
import .ThermalConductivityGroup: ThermalConductivity
import .TemperatureGroup: Temperature
import .RhoCpGroup: RhoCp
import .HeatResidualGroup: HeatResidual

"""
    Arrays <: AbstractArrayCollection

Data structure containing array groups for heat equation solver.

# Fields
- `rhs::`[`RhsHeat`](@ref): Right-hand side arrays for heat equation
- `radiogenic_production::`[`RadiogenicProduction`](@ref): Radiogenic heat 
  production arrays
- `adiabatic_production::`[`AdiabaticProduction`](@ref): Adiabatic heat 
  production arrays
- `thermal_conductivity::`[`ThermalConductivity`](@ref): Thermal conductivity arrays
- `temperature::`[`Temperature`](@ref): Temperature field arrays
- `rhocp::`[`RhoCp`](@ref): Density times heat capacity arrays
- `residuals::`[`HeatResidual`](@ref): Heat equation residual arrays

# Constructor
    Arrays(ynum::Int, xnum::Int)::Arrays

Create a new Arrays collection with the given grid dimensions.

# Arguments
- `ynum::Int`: Number of grid points in y-direction
- `xnum::Int`: Number of grid points in x-direction

"""
mutable struct Arrays <: AbstractArrayCollection
    rhs::RhsHeat
    radiogenic_production::RadiogenicProduction
    adiabatic_production::AdiabaticProduction
    thermal_conductivity::ThermalConductivity
    temperature::Temperature
    rhocp::RhoCp
    residuals::HeatResidual
end

function Arrays(ynum::Int, xnum::Int)::Arrays
    return Arrays(
        RhsHeat(ynum, xnum),
        RadiogenicProduction(ynum, xnum),
        AdiabaticProduction(ynum, xnum),
        ThermalConductivity(ynum, xnum),
        Temperature(ynum, xnum),
        RhoCp(ynum, xnum),
        HeatResidual(ynum, xnum)
    )
end

end # module 