"""
MaterialProperties.jl

Module for managing material property data structures.

# Structs
- `MeltingParameters`: Data structure for storing melting parameters.
- `MaterialColor`: Data structure for storing material color parameters.
- `Density`: Data structure for storing density parameters.
- `DislocationCreep`: Data structure for storing dislocation creep parameters.
- `DiffusionCreep`: Data structure for storing diffusion creep parameters.
- `PeierlsCreep`: Data structure for storing peierls creep parameters.
- `TemperatureDependentViscosity`: Data structure for storing temperature dependent viscosity parameters.
- `Blankenbach89Viscosity`: Data structure for storing blankenbach89 viscosity parameters.
- `ShearModulus`: Data structure for storing shear modulus parameters.
- `Plasticity`: Data structure for storing plasticity parameters.
- `FlowLaw`: Data structure for storing flow law parameters.
- `HeatCapacity`: Data structure for storing heat capacity parameters.
- `ThermalConductivity`: Data structure for storing thermal conductivity parameters.
- `RadiogenicHeatProduction`: Data structure for storing radiogenic heat production parameters.
- `Compaction`: Data structure for storing compaction parameters.

"""
module MaterialProperties

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.SolidusLiquidus.SolidusModels: format_solidus_options
import EarthBox.SolidusLiquidus.LiquidusModels: format_liquidus_options
using EarthBox.Parameters: ParameterFloat, ParameterInt, ParameterStr
using EarthBox.EarthBoxDtypes: MaterialDictType

const PDATA = get_eb_parameters()

function set_value!(
    obj::Any, 
    attribute_name::String, 
    attribute_value::Union{Float64, Int64, String, Nothing}
)::Nothing
    if attribute_value !== nothing
        attribute = getproperty(obj, Symbol(attribute_name))
        if attribute isa ParameterFloat || attribute isa ParameterInt || attribute isa ParameterStr
            attribute.value = attribute_value
        end
    end
    return nothing
end

"""
    MeltingParameters

Data structure for storing melting parameters.

# Fields
- `itype_liquidus::`[`ParameterInt`](@ref): $(PDATA.itype_liquidus.description)
- `stype_liquidus::`[`ParameterStr`](@ref): $(PDATA.stype_liquidus.description)
- `itype_solidus::`[`ParameterInt`](@ref): $(PDATA.itype_solidus.description)
- `stype_solidus::`[`ParameterStr`](@ref): $(PDATA.stype_solidus.description)
- `latent_heat::`[`ParameterFloat`](@ref): $(PDATA.latent_heat.description)

# Returns
- `MeltingParameters`: New MeltingParameters parameter group with initialized values

# Solidus Models
$(format_solidus_options())

# Liquidus Models
$(format_liquidus_options())

"""
mutable struct MeltingParameters
    itype_liquidus::ParameterInt
    stype_liquidus::ParameterStr
    itype_solidus::ParameterInt
    stype_solidus::ParameterStr
    latent_heat::ParameterFloat
end

function set_itype_liquidus!(
    melting_params::MeltingParameters, 
    itype_liquidus::Int64
)::Nothing
    melting_params.itype_liquidus = ParameterInt(
        itype_liquidus, 
        PDATA.itype_liquidus.name, 
        PDATA.itype_liquidus.units,
        PDATA.itype_liquidus.description,
    )
    return nothing
end

function set_itype_solidus!(
    melting_params::MeltingParameters, 
    itype_solidus::Int64
)::Nothing
    melting_params.itype_solidus = ParameterInt(
        itype_solidus, 
        PDATA.itype_solidus.name, 
        PDATA.itype_solidus.units,
        PDATA.itype_solidus.description,
    )
    return nothing
end

function MeltingParameters(input_material_obj::MaterialDictType)::MeltingParameters
    itype_liquidus = ParameterInt(
        input_material_obj[PDATA.itype_liquidus.name],
        PDATA.itype_liquidus.name,
        PDATA.itype_liquidus.units,
        PDATA.itype_liquidus.description,
    )
    stype_liquidus = ParameterStr(
        input_material_obj[PDATA.stype_liquidus.name],
        PDATA.stype_liquidus.name,
        PDATA.stype_liquidus.units,
        PDATA.stype_liquidus.description,
    )
    itype_solidus = ParameterInt(
        input_material_obj[PDATA.itype_solidus.name],
        PDATA.itype_solidus.name,
        PDATA.itype_solidus.units,
        PDATA.itype_solidus.description,
    )
    stype_solidus = ParameterStr(
        input_material_obj[PDATA.stype_solidus.name],
        PDATA.stype_solidus.name,
        PDATA.stype_solidus.units,
        PDATA.stype_solidus.description,
    )
    latent_heat = ParameterFloat(
        input_material_obj[PDATA.latent_heat.name],
        PDATA.latent_heat.name,
        PDATA.latent_heat.units,
        PDATA.latent_heat.description,
    )
    return MeltingParameters(itype_liquidus, stype_liquidus, itype_solidus, stype_solidus, latent_heat)
end

"""
    MaterialColor

Data structure for storing material color parameters.

# Fields
- `red_fraction::`[`ParameterFloat`](@ref): $(PDATA.red_fraction.description)
- `green_fraction::`[`ParameterFloat`](@ref): $(PDATA.green_fraction.description)
- `blue_fraction::`[`ParameterFloat`](@ref): $(PDATA.blue_fraction.description)

# Returns
- `MaterialColor`: New MaterialColor parameter group with initialized values

"""
mutable struct MaterialColor
    red_fraction::ParameterFloat
    green_fraction::ParameterFloat
    blue_fraction::ParameterFloat
end

function MaterialColor(input_material_obj::MaterialDictType)::MaterialColor
    red_fraction = ParameterFloat(
        input_material_obj[PDATA.red_fraction.name],
        PDATA.red_fraction.name,
        PDATA.red_fraction.units,
        PDATA.red_fraction.description,
    )
    green_fraction = ParameterFloat(
        input_material_obj[PDATA.green_fraction.name],
        PDATA.green_fraction.name,
        PDATA.green_fraction.units,
        PDATA.green_fraction.description,
    )
    blue_fraction = ParameterFloat(
        input_material_obj[PDATA.blue_fraction.name],
        PDATA.blue_fraction.name,
        PDATA.blue_fraction.units,
        PDATA.blue_fraction.description,
    )
    return MaterialColor(red_fraction, green_fraction, blue_fraction)
end

"""
    Density

Data structure for storing density parameters.

# Fields
- `standard_density::`[`ParameterFloat`](@ref): $(PDATA.standard_density.description)
- `thermal_expansion::`[`ParameterFloat`](@ref): $(PDATA.thermal_expansion.description)
- `compressibility::`[`ParameterFloat`](@ref): $(PDATA.compressibility.description)
- `melt_density::`[`ParameterFloat`](@ref): $(PDATA.melt_density.description)

# Returns
- `Density`: New Density parameter group with initialized values

"""
mutable struct Density
    standard_density::ParameterFloat
    thermal_expansion::ParameterFloat
    compressibility::ParameterFloat
    melt_density::ParameterFloat
end

function Density(input_material_obj::MaterialDictType)::Density
    standard_density = ParameterFloat(
        input_material_obj[PDATA.standard_density.name],
        PDATA.standard_density.name,
        PDATA.standard_density.units,
        PDATA.standard_density.description,
    )
    thermal_expansion = ParameterFloat(
        input_material_obj[PDATA.thermal_expansion.name],
        PDATA.thermal_expansion.name,
        PDATA.thermal_expansion.units,
        PDATA.thermal_expansion.description,
    )
    compressibility = ParameterFloat(
        input_material_obj[PDATA.compressibility.name],
        PDATA.compressibility.name,
        PDATA.compressibility.units,
        PDATA.compressibility.description,
    )
    melt_density = ParameterFloat(
        input_material_obj[PDATA.melt_density.name],
        PDATA.melt_density.name,
        PDATA.melt_density.units,
        PDATA.melt_density.description,
    )
    return Density(standard_density, thermal_expansion, compressibility, melt_density)
end

"""
    DislocationCreep

Data structure for storing dislocation creep parameters.

# Fields
- `pre_exponential_dc::`[`ParameterFloat`](@ref): $(PDATA.pre_exponential_dc.description)
- `stress_exponent_n_dc::`[`ParameterFloat`](@ref): $(PDATA.stress_exponent_n_dc.description)
- `activation_energy_dc::`[`ParameterFloat`](@ref): $(PDATA.activation_energy_dc.description)
- `activation_volume_dc::`[`ParameterFloat`](@ref): $(PDATA.activation_volume_dc.description)

# Returns
- `DislocationCreep`: New DislocationCreep parameter group with initialized values

"""
mutable struct DislocationCreep
    pre_exponential_dc::ParameterFloat
    stress_exponent_n_dc::ParameterFloat
    activation_energy_dc::ParameterFloat
    activation_volume_dc::ParameterFloat
end

function DislocationCreep(input_material_obj::MaterialDictType)::DislocationCreep
    pre_exponential_dc = ParameterFloat(
        input_material_obj[PDATA.pre_exponential_dc.name],
        PDATA.pre_exponential_dc.name,
        PDATA.pre_exponential_dc.units,
        PDATA.pre_exponential_dc.description,
    )
    stress_exponent_n_dc = ParameterFloat(
        input_material_obj[PDATA.stress_exponent_n_dc.name],
        PDATA.stress_exponent_n_dc.name,
        PDATA.stress_exponent_n_dc.units,
        PDATA.stress_exponent_n_dc.description,
    )
    activation_energy_dc = ParameterFloat(
        input_material_obj[PDATA.activation_energy_dc.name],
        PDATA.activation_energy_dc.name,
        PDATA.activation_energy_dc.units,
        PDATA.activation_energy_dc.description,
    )
    activation_volume_dc = ParameterFloat(
        input_material_obj[PDATA.activation_volume_dc.name],
        PDATA.activation_volume_dc.name,
        PDATA.activation_volume_dc.units,
        PDATA.activation_volume_dc.description,
    )
    return DislocationCreep(pre_exponential_dc, stress_exponent_n_dc, activation_energy_dc, activation_volume_dc)
end

"""
    DiffusionCreep

Data structure for storing diffusion creep parameters.

# Fields
- `pre_exponential_difc::`[`ParameterFloat`](@ref): $(PDATA.pre_exponential_difc.description)
- `activation_energy_difc::`[`ParameterFloat`](@ref): $(PDATA.activation_energy_difc.description)
- `activation_volume_difc::`[`ParameterFloat`](@ref): $(PDATA.activation_volume_difc.description)

# Returns
- `DiffusionCreep`: New DiffusionCreep parameter group with initialized values

"""
mutable struct DiffusionCreep
    pre_exponential_difc::ParameterFloat
    activation_energy_difc::ParameterFloat
    activation_volume_difc::ParameterFloat
end

function DiffusionCreep(input_material_obj::MaterialDictType)::DiffusionCreep
    pre_exponential_difc = ParameterFloat(
        input_material_obj[PDATA.pre_exponential_difc.name],
        PDATA.pre_exponential_difc.name,
        PDATA.pre_exponential_difc.units,
        PDATA.pre_exponential_difc.description,
    )
    activation_energy_difc = ParameterFloat(
        input_material_obj[PDATA.activation_energy_difc.name],
        PDATA.activation_energy_difc.name,
        PDATA.activation_energy_difc.units,
        PDATA.activation_energy_difc.description,
    )
    activation_volume_difc = ParameterFloat(
        input_material_obj[PDATA.activation_volume_difc.name],
        PDATA.activation_volume_difc.name,
        PDATA.activation_volume_difc.units,
        PDATA.activation_volume_difc.description,
    )
    return DiffusionCreep(pre_exponential_difc, activation_energy_difc, activation_volume_difc)
end

"""
    PeierlsCreep

Data structure for storing peierls creep parameters.

# Fields
- `pre_exponential_pei::`[`ParameterFloat`](@ref): $(PDATA.pre_exponential_pei.description)
- `stress_exponent_m1_pei::`[`ParameterFloat`](@ref): $(PDATA.stress_exponent_m1_pei.description)
- `stress_exponent_m2_pei::`[`ParameterFloat`](@ref): $(PDATA.stress_exponent_m2_pei.description)
- `peierls_stress::`[`ParameterFloat`](@ref): $(PDATA.peierls_stress.description)

# Returns
- `PeierlsCreep`: New PeierlsCreep parameter group with initialized values

"""
mutable struct PeierlsCreep
    pre_exponential_pei::ParameterFloat
    stress_exponent_m1_pei::ParameterFloat
    stress_exponent_m2_pei::ParameterFloat
    peierls_stress::ParameterFloat
end

function PeierlsCreep(input_material_obj::MaterialDictType)::PeierlsCreep
    pre_exponential_pei = ParameterFloat(
        input_material_obj[PDATA.pre_exponential_pei.name],
        PDATA.pre_exponential_pei.name,
        PDATA.pre_exponential_pei.units,
        PDATA.pre_exponential_pei.description,
    )
    stress_exponent_m1_pei = ParameterFloat(
        input_material_obj[PDATA.stress_exponent_m1_pei.name],
        PDATA.stress_exponent_m1_pei.name,
        PDATA.stress_exponent_m1_pei.units,
        PDATA.stress_exponent_m1_pei.description,
    )
    stress_exponent_m2_pei = ParameterFloat(
        input_material_obj[PDATA.stress_exponent_m2_pei.name],
        PDATA.stress_exponent_m2_pei.name,
        PDATA.stress_exponent_m2_pei.units,
        PDATA.stress_exponent_m2_pei.description,
    )
    peierls_stress = ParameterFloat(
        input_material_obj[PDATA.peierls_stress.name],
        PDATA.peierls_stress.name,
        PDATA.peierls_stress.units,
        PDATA.peierls_stress.description,
    )
    return PeierlsCreep(pre_exponential_pei, stress_exponent_m1_pei, stress_exponent_m2_pei, peierls_stress)
end

"""
    TemperatureDependentViscosity

Data structure for storing temperature dependent viscosity parameters.

# Fields
- `pre_exponential_td::`[`ParameterFloat`](@ref): $(PDATA.pre_exponential_td.description)
- `activation_energy_td::`[`ParameterFloat`](@ref): $(PDATA.activation_energy_td.description)

# Returns
- `TemperatureDependentViscosity`: New TemperatureDependentViscosity parameter group with initialized values

"""
mutable struct TemperatureDependentViscosity
    pre_exponential_td::ParameterFloat
    activation_energy_td::ParameterFloat
end

function TemperatureDependentViscosity(input_material_obj::MaterialDictType)::TemperatureDependentViscosity
    pre_exponential_td = ParameterFloat(
        input_material_obj[PDATA.pre_exponential_td.name],
        PDATA.pre_exponential_td.name,
        PDATA.pre_exponential_td.units,
        PDATA.pre_exponential_td.description,
    )
    activation_energy_td = ParameterFloat(
        input_material_obj[PDATA.activation_energy_td.name],
        PDATA.activation_energy_td.name,
        PDATA.activation_energy_td.units,
        PDATA.activation_energy_td.description,
    )
    return TemperatureDependentViscosity(pre_exponential_td, activation_energy_td)
end

"""
    Blankenbach89Viscosity

Data structure for storing blankenbach89 viscosity parameters.

# Fields
- `viscosity_ref_blankenbach89::`[`ParameterFloat`](@ref): $(PDATA.viscosity_ref_blankenbach89.description)
- `b_term_blankenbach89::`[`ParameterFloat`](@ref): $(PDATA.b_term_blankenbach89.description)
- `c_term_blankenbach89::`[`ParameterFloat`](@ref): $(PDATA.c_term_blankenbach89.description)

# Returns
- `Blankenbach89Viscosity`: New Blankenbach89Viscosity parameter group with initialized values

"""
mutable struct Blankenbach89Viscosity
    viscosity_ref_blankenbach89::ParameterFloat
    b_term_blankenbach89::ParameterFloat
    c_term_blankenbach89::ParameterFloat
end

function Blankenbach89Viscosity(input_material_obj::MaterialDictType)::Blankenbach89Viscosity
    viscosity_ref_blankenbach89 = ParameterFloat(
        input_material_obj[PDATA.viscosity_ref_blankenbach89.name],
        PDATA.viscosity_ref_blankenbach89.name,
        PDATA.viscosity_ref_blankenbach89.units,
        PDATA.viscosity_ref_blankenbach89.description,
    )
    b_term_blankenbach89 = ParameterFloat(
        input_material_obj[PDATA.b_term_blankenbach89.name],
        PDATA.b_term_blankenbach89.name,
        PDATA.b_term_blankenbach89.units,
        PDATA.b_term_blankenbach89.description,
    )
    c_term_blankenbach89 = ParameterFloat(
        input_material_obj[PDATA.c_term_blankenbach89.name],
        PDATA.c_term_blankenbach89.name,
        PDATA.c_term_blankenbach89.units,
        PDATA.c_term_blankenbach89.description,
    )
    return Blankenbach89Viscosity(viscosity_ref_blankenbach89, b_term_blankenbach89, c_term_blankenbach89)
end

"""
    ShearModulus

Data structure for storing shear modulus parameters.

# Fields
- `shear_modulus::`[`ParameterFloat`](@ref): $(PDATA.shear_modulus.description)

# Returns
- `ShearModulus`: New ShearModulus parameter group with initialized values

"""
mutable struct ShearModulus
    shear_modulus::ParameterFloat
end

function ShearModulus(input_material_obj::MaterialDictType)::ShearModulus
    shear_modulus = ParameterFloat(
        input_material_obj[PDATA.shear_modulus.name],
        PDATA.shear_modulus.name,
        PDATA.shear_modulus.units,
        PDATA.shear_modulus.description,
    )
    return ShearModulus(shear_modulus)
end

"""
    Plasticity

Data structure for storing plasticity parameters.

# Fields
- `cohesion_initial::`[`ParameterFloat`](@ref): $(PDATA.cohesion_initial.description)
- `cohesion_final::`[`ParameterFloat`](@ref): $(PDATA.cohesion_final.description)
- `friction_angle_initial::`[`ParameterFloat`](@ref): $(PDATA.friction_angle_initial.description)
- `friction_angle_final::`[`ParameterFloat`](@ref): $(PDATA.friction_angle_final.description)
- `sine_friction_angle_initial::`[`ParameterFloat`](@ref): $(PDATA.sine_friction_angle_initial.description)
- `sine_friction_angle_final::`[`ParameterFloat`](@ref): $(PDATA.sine_friction_angle_final.description)
- `strain_initial::`[`ParameterFloat`](@ref): $(PDATA.strain_initial.description)
- `strain_final::`[`ParameterFloat`](@ref): $(PDATA.strain_final.description)
- `dilatation_angle::`[`ParameterFloat`](@ref): $(PDATA.dilatation_angle.description)

# Returns
- `Plasticity`: New Plasticity parameter group with initialized values

"""
mutable struct Plasticity
    cohesion_initial::ParameterFloat
    cohesion_final::ParameterFloat
    friction_angle_initial::ParameterFloat
    friction_angle_final::ParameterFloat
    sine_friction_angle_initial::ParameterFloat
    sine_friction_angle_final::ParameterFloat
    strain_initial::ParameterFloat
    strain_final::ParameterFloat
    dilatation_angle::ParameterFloat
end

function Plasticity(input_material_obj::MaterialDictType)::Plasticity
    friction_angle_initial = ParameterFloat(
        input_material_obj[PDATA.friction_angle_initial.name],
        PDATA.friction_angle_initial.name,
        PDATA.friction_angle_initial.units,
        PDATA.friction_angle_initial.description,
    )
    friction_angle_final = ParameterFloat(
        input_material_obj[PDATA.friction_angle_final.name],
        PDATA.friction_angle_final.name,
        PDATA.friction_angle_final.units,
        PDATA.friction_angle_final.description,
    )
    
    # Calculate sine values
    sine_friction_angle_initial = ParameterFloat(
        sin(friction_angle_initial.value * π / 180.0),
        PDATA.sine_friction_angle_initial.name,
        PDATA.sine_friction_angle_initial.units,
        PDATA.sine_friction_angle_initial.description,
    )
    sine_friction_angle_final = ParameterFloat(
        sin(friction_angle_final.value * π / 180.0),
        PDATA.sine_friction_angle_final.name,
        PDATA.sine_friction_angle_final.units,
        PDATA.sine_friction_angle_final.description,
    )
    return Plasticity(
        ParameterFloat(
            input_material_obj[PDATA.cohesion_initial.name],
            PDATA.cohesion_initial.name,
            PDATA.cohesion_initial.units,
            PDATA.cohesion_initial.description,
        ),
        ParameterFloat(
            input_material_obj[PDATA.cohesion_final.name],
            PDATA.cohesion_final.name,
            PDATA.cohesion_final.units,
            PDATA.cohesion_final.description,
        ),
        friction_angle_initial,
        friction_angle_final,
        sine_friction_angle_initial,
        sine_friction_angle_final,
        ParameterFloat(
            input_material_obj[PDATA.strain_initial.name],
            PDATA.strain_initial.name,
            PDATA.strain_initial.units,
            PDATA.strain_initial.description,
        ),
        ParameterFloat(
            input_material_obj[PDATA.strain_final.name],
            PDATA.strain_final.name,
            PDATA.strain_final.units,
            PDATA.strain_final.description,
        ),
        ParameterFloat(
            input_material_obj[PDATA.dilatation_angle.name],
            PDATA.dilatation_angle.name,
            PDATA.dilatation_angle.units,
            PDATA.dilatation_angle.description,
        )
    )
end

"""
    FlowLaw

Data structure for storing flow law parameters.

# Fields
- `flow_type::`[`ParameterInt`](@ref): $(PDATA.flow_type.description)
- `flow_stype::`[`ParameterStr`](@ref): $(PDATA.flow_stype.description)
- `viscosity_iso::`[`ParameterFloat`](@ref): $(PDATA.viscosity_iso.description)
- `dislocation_creep::`[`DislocationCreep`](@ref): DislocationCreep parameter group
- `peierls_creep::`[`PeierlsCreep`](@ref): PeierlsCreep parameter group
- `diffusion_creep::`[`DiffusionCreep`](@ref): DiffusionCreep parameter group
- `temperature_dependent_viscosity::`[`TemperatureDependentViscosity`](@ref): TemperatureDependentViscosity parameter group
- `blankenbach89_viscosity::`[`Blankenbach89Viscosity`](@ref): Blankenbach89Viscosity parameter group
- `plasticity::`[`Plasticity`](@ref): Plasticity parameter group
- `shear_modulus::`[`ShearModulus`](@ref): ShearModulus parameter group

# Returns
- `FlowLaw`: New FlowLaw parameter group with initialized values

"""
mutable struct FlowLaw
    flow_type::ParameterInt
    flow_stype::ParameterStr
    viscosity_iso::ParameterFloat
    dislocation_creep::DislocationCreep
    peierls_creep::PeierlsCreep
    diffusion_creep::DiffusionCreep
    temperature_dependent_viscosity::TemperatureDependentViscosity
    blankenbach89_viscosity::Blankenbach89Viscosity
    plasticity::Plasticity
    shear_modulus::ShearModulus
end

function FlowLaw(input_material_obj::MaterialDictType)::FlowLaw
    flow_stype = ParameterStr(
        input_material_obj[PDATA.flow_stype.name],
        PDATA.flow_stype.name,
        PDATA.flow_stype.units,
        PDATA.flow_stype.description,
    )
    flow_type = ParameterInt(
        input_material_obj[PDATA.flow_type.name],
        PDATA.flow_type.name,
        PDATA.flow_type.units,
        PDATA.flow_type.description,
    )
    viscosity_iso = ParameterFloat(
        input_material_obj[PDATA.viscosity_iso.name],
        PDATA.viscosity_iso.name,
        PDATA.viscosity_iso.units,
        PDATA.viscosity_iso.description,
    )
    
    return FlowLaw(
        flow_type,
        flow_stype,
        viscosity_iso,
        DislocationCreep(input_material_obj),
        PeierlsCreep(input_material_obj),
        DiffusionCreep(input_material_obj),
        TemperatureDependentViscosity(input_material_obj),
        Blankenbach89Viscosity(input_material_obj),
        Plasticity(input_material_obj),
        ShearModulus(input_material_obj)
    )
end

function set_flow_type!(flow_law::FlowLaw, flow_type::Int64)::Nothing
    flow_law.flow_type = ParameterInt(
        flow_type,
        PDATA.flow_type.name,
        PDATA.flow_type.units,
        PDATA.flow_type.description,
    )
    return nothing
end

function get_flow_type_descriptions()::Dict{String, String}
    flow_dict = Dict{String, String}()
    flow_dict["-1"] = "Viscosity as a function of x"
    flow_dict["0"] = "Isoviscous"
    flow_dict["1"] = "Composite rheology"
    flow_dict["4"] = "Couette flow benchmark"
    flow_dict["5"] = "Convection in a box (Blankenbach et al., 1989)"
    return flow_dict
end

"""
    HeatCapacity

Data structure for storing heat capacity parameters.

# Fields
- `heat_capacity::`[`ParameterFloat`](@ref): $(PDATA.heat_capacity.description)

# Returns
- `HeatCapacity`: New HeatCapacity parameter group with initialized values

"""
mutable struct HeatCapacity
    heat_capacity::ParameterFloat
end

function HeatCapacity(input_material_obj::MaterialDictType)::HeatCapacity
    heat_capacity = ParameterFloat(
        input_material_obj[PDATA.heat_capacity.name],
        PDATA.heat_capacity.name,
        PDATA.heat_capacity.units,
        PDATA.heat_capacity.description,
    )
    return HeatCapacity(heat_capacity)
end

"""
    ThermalConductivity

Data structure for storing thermal conductivity parameters.

# Fields
- `thermal_conductivity_ref::`[`ParameterFloat`](@ref): $(PDATA.thermal_conductivity_ref.description)
- `thermal_conductivity_a::`[`ParameterFloat`](@ref): $(PDATA.thermal_conductivity_a.description)

# Returns
- `ThermalConductivity`: New ThermalConductivity parameter group with initialized values

"""
mutable struct ThermalConductivity
    thermal_conductivity_ref::ParameterFloat
    thermal_conductivity_a::ParameterFloat
end

function ThermalConductivity(input_material_obj::MaterialDictType)::ThermalConductivity
    thermal_conductivity_ref = ParameterFloat(
        input_material_obj[PDATA.thermal_conductivity_ref.name],
        PDATA.thermal_conductivity_ref.name,
        PDATA.thermal_conductivity_ref.units,
        PDATA.thermal_conductivity_ref.description,
    )
    thermal_conductivity_a = ParameterFloat(
        input_material_obj[PDATA.thermal_conductivity_a.name],
        PDATA.thermal_conductivity_a.name,
        PDATA.thermal_conductivity_a.units,
        PDATA.thermal_conductivity_a.description,
    )
    return ThermalConductivity(thermal_conductivity_ref, thermal_conductivity_a)
end

"""
    RadiogenicHeatProduction

Data structure for storing radiogenic heat production parameters.

# Fields
- `radiogenic_heat_production::`[`ParameterFloat`](@ref): $(PDATA.radiogenic_heat_production.description)

# Returns
- `RadiogenicHeatProduction`: New RadiogenicHeatProduction parameter group with initialized values

"""
mutable struct RadiogenicHeatProduction
    radiogenic_heat_production::ParameterFloat
end

function RadiogenicHeatProduction(input_material_obj::MaterialDictType)::RadiogenicHeatProduction
    radiogenic_heat_production = ParameterFloat(
        input_material_obj[PDATA.radiogenic_heat_production.name],
        PDATA.radiogenic_heat_production.name,
        PDATA.radiogenic_heat_production.units,
        PDATA.radiogenic_heat_production.description,
    )
    return RadiogenicHeatProduction(radiogenic_heat_production)
end

"""
    Compaction

Data structure for storing compaction parameters.

# Fields
- `porosity_initial::`[`ParameterFloat`](@ref): $(PDATA.porosity_initial.description)
- `porosity_decay_depth::`[`ParameterFloat`](@ref): $(PDATA.porosity_decay_depth.description)

# Returns
- `Compaction`: New Compaction parameter group with initialized values

"""
mutable struct Compaction
    porosity_initial::ParameterFloat
    porosity_decay_depth::ParameterFloat
end

function Compaction(input_material_obj::MaterialDictType)::Compaction
    porosity_initial = ParameterFloat(
        input_material_obj[PDATA.porosity_initial.name],
        PDATA.porosity_initial.name,
        PDATA.porosity_initial.units,
        PDATA.porosity_initial.description,
    )
    porosity_decay_depth = ParameterFloat(
        input_material_obj[PDATA.porosity_decay_depth.name],
        PDATA.porosity_decay_depth.name,
        PDATA.porosity_decay_depth.units,
        PDATA.porosity_decay_depth.description,
    )
    return Compaction(porosity_initial, porosity_decay_depth)
end

"""
    MaterialPropertiesState

Data structure for storing material properties state.

# Fields
- `melting_params::`[`MeltingParameters`](@ref): Melting parameters
- `material_color::`[`MaterialColor`](@ref): Material color
- `density::`[`Density`](@ref): Density
- `dislocation_creep::`[`DislocationCreep`](@ref): Dislocation creep
- `diffusion_creep::`[`DiffusionCreep`](@ref): Diffusion creep
- `peierls_creep::`[`PeierlsCreep`](@ref): Peierls creep
- `temperature_dependent_viscosity::`[`TemperatureDependentViscosity`](@ref): Temperature dependent viscosity
- `blankenbach89_viscosity::`[`Blankenbach89Viscosity`](@ref): Blankenbach89 viscosity
- `shear_modulus::`[`ShearModulus`](@ref): Shear modulus
- `plasticity::`[`Plasticity`](@ref): Plasticity
- `flow_law::`[`FlowLaw`](@ref): Flow law
- `heat_capacity::`[`HeatCapacity`](@ref): Heat capacity
- `thermal_conductivity::`[`ThermalConductivity`](@ref): Thermal conductivity
- `radiogenic_heat_production::`[`RadiogenicHeatProduction`](@ref): Radiogenic heat production
- `compaction::`[`Compaction`](@ref): Compaction

# Returns
- `MaterialPropertiesState`: New MaterialPropertiesState parameter group with initialized values

---
# List of Material Properties For Material Collection Files

For a given material in a material collection library file enter the name of the properties followed 
by a list containing value, units and description:

```
MaterialName:
    property_name_1: [value, units, description]
    property_name_2: [value, units, description]
```
Use the standard units defined below. A future version of EarthBox will support material property 
unit conversion.

The following properties can be entered for a given material:

[`MeltingParameters`](@ref)
- `itype_liquidus`: $(PDATA.itype_liquidus.description)
- `stype_liquidus`: $(PDATA.stype_liquidus.description)
- `itype_solidus`: $(PDATA.itype_solidus.description)
- `stype_solidus`: $(PDATA.stype_solidus.description)
- `latent_heat`: $(PDATA.latent_heat.description)
[`Density`](@ref)
- `standard_density`: $(PDATA.standard_density.description)
- `thermal_expansion`: $(PDATA.thermal_expansion.description)
- `compressibility`: $(PDATA.compressibility.description)
- `melt_density`: $(PDATA.melt_density.description)
[`DislocationCreep`](@ref)
- `pre_exponential_dc`: $(PDATA.pre_exponential_dc.description)
- `stress_exponent_n_dc`: $(PDATA.stress_exponent_n_dc.description)
- `activation_energy_dc`: $(PDATA.activation_energy_dc.description)
- `activation_volume_dc`: $(PDATA.activation_volume_dc.description)
[`DiffusionCreep`](@ref)
- `pre_exponential_difc`: $(PDATA.pre_exponential_difc.description)
- `activation_energy_difc`: $(PDATA.activation_energy_difc.description)
- `activation_volume_difc`: $(PDATA.activation_volume_difc.description)
[`PeierlsCreep`](@ref)
- `pre_exponential_pei`: $(PDATA.pre_exponential_pei.description)
- `stress_exponent_m1_pei`: $(PDATA.stress_exponent_m1_pei.description)
- `stress_exponent_m2_pei`: $(PDATA.stress_exponent_m2_pei.description)
- `peierls_stress`: $(PDATA.peierls_stress.description)
[`TemperatureDependentViscosity`](@ref)
- `pre_exponential_td`: $(PDATA.pre_exponential_td.description)
- `activation_energy_td`: $(PDATA.activation_energy_td.description)
[`Blankenbach89Viscosity`](@ref)
- `viscosity_ref_blankenbach89`: $(PDATA.viscosity_ref_blankenbach89.description)
- `b_term_blankenbach89`: $(PDATA.b_term_blankenbach89.description)
- `c_term_blankenbach89`: $(PDATA.c_term_blankenbach89.description)
[`ShearModulus`](@ref)
- `shear_modulus`: $(PDATA.shear_modulus.description)
[`Plasticity`](@ref)
- `cohesion_initial`: $(PDATA.cohesion_initial.description)
- `cohesion_final`: $(PDATA.cohesion_final.description)
- `friction_angle_initial`: $(PDATA.friction_angle_initial.description)
- `friction_angle_final`: $(PDATA.friction_angle_final.description)
- `strain_initial`: $(PDATA.strain_initial.description)
- `strain_final`: $(PDATA.strain_final.description)
- `dilatation_angle`: $(PDATA.dilatation_angle.description)
[`FlowLaw`](@ref)
- `flow_type`: $(PDATA.flow_type.description)
- `flow_stype`: $(PDATA.flow_stype.description)
[`HeatCapacity`](@ref)
- `heat_capacity`: $(PDATA.heat_capacity.description)
[`ThermalConductivity`](@ref)
- `thermal_conductivity_ref`: $(PDATA.thermal_conductivity_ref.description)
- `thermal_conductivity_a`: $(PDATA.thermal_conductivity_a.description)
[`RadiogenicHeatProduction`](@ref)
- `radiogenic_heat_production`: $(PDATA.radiogenic_heat_production.description)
[`Compaction`](@ref)
- `porosity_initial`: $(PDATA.porosity_initial.description)
- `porosity_decay_depth`: $(PDATA.porosity_decay_depth.description)

"""
mutable struct MaterialPropertiesState
    melting_params::MeltingParameters
    material_color::MaterialColor
    density::Density
    dislocation_creep::DislocationCreep
    diffusion_creep::DiffusionCreep
    peierls_creep::PeierlsCreep
    temperature_dependent_viscosity::TemperatureDependentViscosity
    blankenbach89_viscosity::Blankenbach89Viscosity
    shear_modulus::ShearModulus
    plasticity::Plasticity
    flow_law::FlowLaw
    heat_capacity::HeatCapacity
    thermal_conductivity::ThermalConductivity
    radiogenic_heat_production::RadiogenicHeatProduction
    compaction::Compaction
end

function MaterialPropertiesState(input_material_obj::MaterialDictType)::MaterialPropertiesState
    melting_params = MeltingParameters(input_material_obj)
    material_color = MaterialColor(input_material_obj)
    density = Density(input_material_obj)
    dislocation_creep = DislocationCreep(input_material_obj)
    diffusion_creep = DiffusionCreep(input_material_obj)
    peierls_creep = PeierlsCreep(input_material_obj)
    temperature_dependent_viscosity = TemperatureDependentViscosity(input_material_obj)
    blankenbach89_viscosity = Blankenbach89Viscosity(input_material_obj)
    shear_modulus = ShearModulus(input_material_obj)
    plasticity = Plasticity(input_material_obj)
    flow_law = FlowLaw(input_material_obj)
    heat_capacity = HeatCapacity(input_material_obj)
    thermal_conductivity = ThermalConductivity(input_material_obj)
    radiogenic_heat_production = RadiogenicHeatProduction(input_material_obj)
    compaction = Compaction(input_material_obj)
    return MaterialPropertiesState(
        melting_params, material_color, density, dislocation_creep, diffusion_creep, 
        peierls_creep, temperature_dependent_viscosity, blankenbach89_viscosity, 
        shear_modulus, plasticity, flow_law, heat_capacity, thermal_conductivity, 
        radiogenic_heat_production, compaction
        )
end


end # module MaterialProperties 