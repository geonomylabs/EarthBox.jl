module ArrayCollection

import EarthBox.PrintFuncs: print_info
import EarthBox.EarthBoxDtypes
import EarthBox.Arrays.ArrayTypes.MaterialArrayFloat1D: MaterialArrayFloat1DState
import EarthBox.Arrays.ArrayTypes.MaterialArrayFloat2D: MaterialArrayFloat2DState
import EarthBox.Arrays.ArrayTypes.MaterialArrayInt1D: MaterialArrayInt1DState
import EarthBox.Arrays.ArrayTypes.MaterialArrayInt2D: MaterialArrayInt2DState
import EarthBox.SolidusLiquidus: SolidusModels
import EarthBox.SolidusLiquidus: LiquidusModels
import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: get_id
import EarthBox.OptionTools: get_name

const ROOT_NAME = "model.materials.arrays"

"""
    Arrays

Collection of material property arrays used in EarthBox.

# Fields
- `mat_rho::`[`MaterialArrayFloat2DState`](@ref): (nmats, 4) : Density-related properties
  - Column 1: Standard density (kg/m^3)
  - Column 2: Thermal expansion coefficient (1/K)
  - Column 3: Compressibility (1/Pa)
  - Column 4: Melt density (kg/m^3)
- `mat_flow_type::`[`MaterialArrayInt1DState`](@ref): (nmats) : Flow law type integer ID for each material
- `mat_flow::`[`MaterialArrayFloat2DState`](@ref): (nmats, 17) : Flow law parameters
  - Column 1: Viscosity for isoviscous models (Pa.s)
  - Column 2: Pre-exponential factor for power law (1/s/MPa^n)
  - Column 3: Power law stress exponent n
  - Column 4: Activation energy for power law (kJ/mol)
  - Column 5: Activation volume for power law (J/MPa/mol)
  - Column 6: Pre-exponential term for Peierls creep (s^-1*MPa^-2)
  - Column 7: Peierls creep stress exponent 1
  - Column 8: Peierls creep stress exponent 2
  - Column 9: Peierls stress (MPa)
  - Column 10: Pre-exponential terms for diffusion creep (1/s/MPa^n)
  - Column 11: Activation energy for diffusion creep (kJ/mol)
  - Column 12: Activation volume for diffusion creep (J/MPa/mol)
  - Column 13: Pre-exponential term for temperature-dependent viscosity (Pa.s)
  - Column 14: Activation energy for temperature-dependent viscosity (kJ/mol)
  - Column 15: Reference viscosity for flow law from Blankenbach (1989) (Pa.s)
  - Column 16: b term for flow law from Blankenbach (1989)
  - Column 17: c term for flow law from Blankenbach (1989)
- `mat_mu::`[`MaterialArrayFloat1DState`](@ref): (nmats) : Shear modulus (Pa)
  - Column 1: Shear modulus (Pa)
- `mat_plastic::`[`MaterialArrayFloat2DState`](@ref): (nmats, 7) : Plastic deformation parameters
  - Column 1: Cohesion (C0) at initial strain GAM0 (Pa)
  - Column 2: Cohesion (C1) at final strain GAM1 (Pa)
  - Column 3: Sine of friction angle at strain GAM0
  - Column 3: Sine of friction angle at strain GAM1
  - Column 5: Initial strain GAM0
  - Column 6: Final strain GAM1
  - Column 7: Dilatation angle
- `mat_cp::`[`MaterialArrayFloat1DState`](@ref): (nmats) : Heat capacity (J/kg/K)
- `mat_kt::`[`MaterialArrayFloat2DState`](@ref): (nmats, 2) : Thermal conductivity parameters
  - Column 1: Reference thermal conductivity k0 (W/m/k)
  - Column 2: Thermal conductivity term a (W/m)
- `mat_hr::`[`MaterialArrayFloat1DState`](@ref): (nmats) : Radiogenic heat production rate (W/m^3)
- `mat_melting_itypes::`[`MaterialArrayInt2DState`](@ref): (nmats, 2) : Melting model type IDs
  - Column 1: Integer type ID for solidus
  - Column 2: Integer type ID for liquidus
- `mat_melting::`[`MaterialArrayFloat2DState`](@ref): (nmats, 2) : Melting model parameters
  - Column 1: Latent heat (J/kg)
  - Column 2: Shift in solidus and liquidus temperature to account for fractional crystallization above the Moho at spreading ridges. (K)
- `mat_compaction::`[`MaterialArrayFloat2DState`](@ref): Compaction parameters
  - Column 1: Initial porosity at sediment-water interface (fraction)
  - Column 2: Porosity decay depth (meters)

# Nested Dot Access
- `mat_rho = $(ROOT_NAME).mat_rho.array`
- `mat_flow_type = $(ROOT_NAME).mat_flow_type.array`
- `mat_flow = $(ROOT_NAME).mat_flow.array`
- `mat_mu = $(ROOT_NAME).mat_mu.array`
- `mat_plastic = $(ROOT_NAME).mat_plastic.array`
- `mat_cp = $(ROOT_NAME).mat_cp.array`
- `mat_kt = $(ROOT_NAME).mat_kt.array`
- `mat_hr = $(ROOT_NAME).mat_hr.array`
- `mat_melting_itypes = $(ROOT_NAME).mat_melting_itypes.array`
- `mat_melting = $(ROOT_NAME).mat_melting.array`
- `mat_compaction = $(ROOT_NAME).mat_compaction.array`

# Constructor
    Arrays(nmats::Int)

# Arguments
- `nmats::Int`: Number of materials

# Returns
- `Arrays`: New Arrays collection with initialized arrays
"""
mutable struct Arrays
    mat_rho::MaterialArrayFloat2DState
    mat_flow_type::MaterialArrayInt1DState
    mat_flow::MaterialArrayFloat2DState
    mat_mu::MaterialArrayFloat1DState
    mat_plastic::MaterialArrayFloat2DState
    mat_cp::MaterialArrayFloat1DState
    mat_kt::MaterialArrayFloat2DState
    mat_hr::MaterialArrayFloat1DState
    mat_melting_itypes::MaterialArrayInt2DState
    mat_melting::MaterialArrayFloat2DState
    mat_compaction::MaterialArrayFloat2DState
end

function Arrays(nmats::Int)
    return Arrays(
        MaterialArrayFloat2DState(
            zeros(nmats, 4),
            "mat_rho",
            [
                "1: kg/m^3",
                "2: 1/K",
                "3: 1/Pa",
                "4: kg/m^3"
            ],
            [
                "1: Standard density",
                "2: Thermal expansion",
                "3: Compressibility",
                "4: Melt density"
            ]
        ),
        MaterialArrayInt1DState(
            fill(Int64(-9999), nmats),
            "mat_flow_type",
            [
                "None"
            ],
            [
                ("Flow law type integer ID. See get_option() method in "
                *"rheology.flow_viscosity for details on available options.")
            ]
        ),
        MaterialArrayFloat2DState(
            zeros(nmats, 17),
            "mat_flow",
            [
                "1: Pa.s",
                "2: 1/s/MPa^n",
                "3: None",
                "4: kJ/mol",
                "5: J/MPa/mol",
                "6: s^-1*MPa^-2",
                "7: None",
                "8: None",
                "9: MPa",
                "10: 1/s/MPa^n",
                "11: kJ/mol",
                "12: J/MPa/mol",
                "13: Pa.s",
                "14: kJ/mol",
                "15: Pa.s",
                "16: None",
                "17: None"
            ],
            [
                "1: Viscosity for isoviscous models",
                "2: Pre-exponential factor for power law",
                "3: Power law stress exponent n",
                "4: Activation energy for power law",
                "5: Activation volume for power law",
                "6: Pre-exponential term for Peierls creep",
                "7: Peierls creep stress exponent 1",
                "8: Peierls creep stress exponent 2",
                "9: Peierls stress",
                "10: Pre-exponential terms for diffusion creep",
                "11: Activation energy for diffusion creep",
                "12: Activation volume for diffusion creep",
                "13: Pre-exponential term for temperature-dependent viscosity",
                "14: Activation energy for temperature-dependent viscosity",
                "15: Reference viscosity for flow law from Blankenbach (1989).",
                "16: b term for flow law from Blankenbach (1989).",
                "17: c term for flow law from Blankenbach (1989)."
            ]
        ),
        MaterialArrayFloat1DState(
            zeros(nmats),
            "mat_mu",
            ["1: Pa"],
            ["1: Shear modulus"]
        ),
        MaterialArrayFloat2DState(
            zeros(nmats, 7),
            "mat_plastic",
            [
                "1: Pa",
                "2: Pa",
                "3: None",
                "4: None",
                "5: None",
                "6: None",
                "7: Degrees"
            ],
            [
                "1: Cohesion (C0) at initial strain GAM0",
                "2: Cohesion (C1) at final strain GAM1",
                "3: Sine of friction angle at strain GAM0",
                "4: Sine of friction angle at strain GAM1",
                "5: Initial strain GAM0",
                "6: Final strain GAM1",
                "7: Dilatation angle"
            ]
        ),
        MaterialArrayFloat1DState(
            zeros(nmats),
            "mat_cp",
            ["1: J/kg/K"],
            ["1: Heat capacity."]
        ),
        MaterialArrayFloat2DState(
            zeros(nmats, 2),
            "mat_kt",
            [
                "1: W/m/k",
                "2: W/m"
            ],
            [
                "1: Reference thermal conductivity k0",
                "2: Thermal conductivity term a"
            ]
        ),
        MaterialArrayFloat1DState(
            zeros(nmats),
            "mat_hr",
            ["1: W/m^3"],
            ["2: Radiogenic heat production"]
        ),
        MaterialArrayInt2DState(
            zeros(Int64, nmats, 2),
            "mat_melting_itypes",
            [
                "1: None",
                "2: None"
            ],
            [
                "1: Integer type ID for solidus",
                "2: Integer type ID for liquidus"
            ]
        ),
        MaterialArrayFloat2DState(
            zeros(nmats, 2),
            "mat_melting",
            [
                "1: J/kg",
                "2: K"
            ],
            [
                "1: Latent heat",
                "2: Shift in solidus and liquidus temperature to account for "
                *"fractional crystallization above the Moho at spreading ridges."
            ]
        ),
        MaterialArrayFloat2DState(
            zeros(nmats, 2),
            "mat_compaction",
            [
                "1: fraction",
                "2: meters"
            ],
            [
                "1: Initial porosity at sediment-water interface",
                "2: Porosity decay depth"
            ]
        )
    )
end

function set_cohesion_at_minimum_strain(
    state::Arrays,
    material_id::Int16,
    cohesion::Union{Float64, Nothing}=nothing,
    print_message::Bool=true
)
    if isnothing(cohesion)
        return
    end
    old_value = state.mat_plastic.array[material_id, 1]
    new_value = cohesion
    state.mat_plastic.array[material_id, 1] = new_value
    print_override(
        print_message, old_value, new_value, material_id, "cohesion initial")
end

function set_cohesion_at_maximum_strain(
    state::Arrays,
    material_id::Int16,
    cohesion::Union{Float64, Nothing}=nothing,
    print_message::Bool=true
)
    if isnothing(cohesion)
        return
    end
    old_value = state.mat_plastic.array[material_id, 2]
    new_value = cohesion
    state.mat_plastic.array[material_id, 2] = new_value
    print_override(
        print_message, old_value, new_value, material_id, "cohesion final"
        )
end

function set_friction_angle_initial(
    state::Arrays,
    material_id::Int16,
    friction_angle_initial::Union{Float64, Nothing}=nothing,
    print_message::Bool=true
)
    if isnothing(friction_angle_initial)
        return
    end
    old_value = state.mat_plastic.array[material_id, 3]
    new_value = sin(deg2rad(friction_angle_initial))
    state.mat_plastic.array[material_id, 3] = new_value
    print_override(
        print_message, old_value, new_value, material_id, 
        "friction coefficient initial"
        )
end

function set_friction_angle_final(
    state::Arrays,
    material_id::Int16,
    friction_angle_final::Union{Float64, Nothing}=nothing,
    print_message::Bool=true
)
    if isnothing(friction_angle_final)
        return
    end
    old_value = state.mat_plastic.array[material_id, 4]
    new_value = sin(deg2rad(friction_angle_final))
    state.mat_plastic.array[material_id, 4] = new_value
    print_override(
        print_message, old_value, new_value, material_id, 
        "friction coefficient final"
        )
end

function set_initial_strain(
    state::Arrays,
    material_id::Int,
    initial_strain::Union{Float64, Nothing}=nothing,
    print_message::Bool=true
)
    if isnothing(initial_strain)
        return
    end
    old_value = state.mat_plastic.array[material_id, 5]
    new_value = initial_strain
    state.mat_plastic.array[material_id, 5] = new_value
    print_override(
        print_message, old_value, new_value, material_id, "initial strain")
end

function set_final_strain(
    state::Arrays,
    material_id::Int16,
    strain_final::Union{Float64, Nothing}=nothing,
    print_message::Bool=true
)
    if isnothing(strain_final)
        return
    end
    old_value = state.mat_plastic.array[material_id, 6]
    new_value = strain_final
    state.mat_plastic.array[material_id, 6] = new_value
    print_override(
        print_message, old_value, new_value, material_id, "final strain")
end

function scale_pre_exponential_factor_dislocation_creep(
    state::Arrays,
    material_id::Int16,
    scale_factor::Union{Float64, Nothing}=nothing,
    print_message::Bool=true
)
    if isnothing(scale_factor)
        return
    end
    old_value = state.mat_flow.array[material_id, 2]
    new_value = old_value * scale_factor
    state.mat_flow.array[material_id, 2] = new_value
    print_override(
        print_message, old_value, new_value, material_id, 
        "pre-exponential factor dislocation creep"
        )
end

function scale_pre_exponential_factor_diffusion_creep(
    state::Arrays,
    material_id::Int16,
    scale_factor::Union{Float64, Nothing}=nothing,
    print_message::Bool=true
)
    if isnothing(scale_factor)
        return
    end
    old_value = state.mat_flow.array[material_id, 10]
    new_value = old_value * scale_factor
    state.mat_flow.array[material_id, 10] = new_value
    print_override(
        print_message, old_value, new_value, material_id, 
        "pre-exponential factor diffusion creep"
        )
end

function set_radiogenic_heat_production(
    state::Arrays,
    material_id::Int16,
    heat_production::Union{Float64, Nothing}=nothing,
    print_message::Bool=true
)
    if isnothing(heat_production)
        return
    end
    old_value = state.mat_hr.array[material_id]
    new_value = heat_production
    state.mat_hr.array[material_id] = new_value
    print_override(
        print_message, old_value, new_value, material_id, 
        "radiogenic heat production"
        )
end

function set_latent_heat(
    state::Arrays,
    material_id::Int16,
    latent_heat::Union{Float64, Nothing}=nothing,
    print_message::Bool=true
)
    if isnothing(latent_heat)
        return
    end
    old_value = state.mat_melting.array[material_id, 1]
    new_value = latent_heat
    state.mat_melting.array[material_id, 1] = new_value
    print_override(print_message, old_value, new_value, material_id, "latent heat")
end

function set_solidus_model(
    state::Arrays,
    material_id::Int16,
    solidus_model_new::Union{String, Nothing}=nothing,
    print_message::Bool=true
)::Nothing
    if isnothing(solidus_model_new)
        return
    end
    options_dict = SolidusModels.get_options()
    
    itype_solidus_old = state.mat_melting_itypes.array[material_id, 1]
    solidus_model_old = get_name(options_dict, itype_solidus_old)

    itype_solidus_new = get_id(options_dict, solidus_model_new)
    state.mat_melting_itypes.array[material_id, 1] = itype_solidus_new
    
    print_override(print_message, solidus_model_old, solidus_model_new, material_id, "solidus model")
    return nothing
end

function set_liquidus_model(
    state::Arrays,
    material_id::Int16,
    liquidus_model_new::Union{String, Nothing}=nothing,
    print_message::Bool=true
)::Nothing
    if isnothing(liquidus_model_new)
        return
    end
    options_dict = LiquidusModels.get_options()
    
    itype_liquidus_old = state.mat_melting_itypes.array[material_id, 2]
    liquidus_model_old = get_name(options_dict, itype_liquidus_old)
    
    itype_liquidus_new = get_id(options_dict, liquidus_model_new)
    state.mat_melting_itypes.array[material_id, 2] = itype_liquidus_new
    
    print_override(print_message, liquidus_model_old, liquidus_model_new, material_id, "liquidus model")
    return nothing
end

function print_override(
    print_message::Bool,
    old_value::Union{Float64, Int64, String},
    new_value::Union{Float64, Int64, String},
    material_id::Int16,
    property_name::String
)
    if print_message && old_value != new_value
        print_info(
            ">> Overriding $property_name for material $material_id from $old_value to $new_value."
            )
    end
end

end # module 