module MaterialContainer

import EarthBox.PrintFuncs: print_info
import EarthBox.Parameters: ParameterInt
import EarthBox.Parameters: ParameterStr
import EarthBox.EarthBoxDtypes: MaterialDictType, MaterialsDictType
import ..MaterialProperties: Density
import ..MaterialProperties: FlowLaw
import ..MaterialProperties: ShearModulus
import ..MaterialProperties: Plasticity
import ..MaterialProperties: HeatCapacity
import ..MaterialProperties: ThermalConductivity
import ..MaterialProperties: RadiogenicHeatProduction
import ..MaterialProperties: MaterialColor
import ..MaterialProperties: MeltingParameters
import ..MaterialProperties: Compaction
import ...Registry: get_types_registry
import ...Registry: get_domains_registry
import ...Registry: MaterialDomainsRegistry

mutable struct Material
    matid::ParameterInt
    mat_name::ParameterStr
    mat_type::ParameterStr
    mat_domain::ParameterStr
    description::String
    density::Density
    flow_law::FlowLaw
    shear_modulus::ShearModulus
    plasticity::Plasticity
    heat_capacity::HeatCapacity
    thermal_conductivity::ThermalConductivity
    radiogenic_heatproduction::RadiogenicHeatProduction
    melting_parameters::MeltingParameters
    compaction::Compaction
    rgb::MaterialColor
end

function Material(material_dict::MaterialDictType)::Material
    if isa(material_dict["matid"], String)
        matid = parse(Int16, material_dict["matid"])
    else
        matid = material_dict["matid"]  # Keep as Int16
    end
    matid = ParameterInt(
        convert(Int64, matid),
        "matid",
        "None",
        "Integer ID of material."
    )
    mat_name = ParameterStr(
        material_dict["mat_name"],
        "mat_name",
        "None",
        "Name of material from library."
    )
    mat_type = ParameterStr(
        get_type_name(material_dict),
        "mat_type",
        "None",
        "Genetic type of material."
    )
    mat_domain = ParameterStr(
        get_domain_name(material_dict),
        "mat_domain",
        "None",
        "Domain of material."
    )
    
    return Material(
        matid,
        mat_name,
        mat_type,
        mat_domain,
        "",
        Density(material_dict),
        FlowLaw(material_dict),
        ShearModulus(material_dict),
        Plasticity(material_dict),
        HeatCapacity(material_dict),
        ThermalConductivity(material_dict),
        RadiogenicHeatProduction(material_dict),
        MeltingParameters(material_dict),
        Compaction(material_dict),
        MaterialColor(material_dict)
    )
end

function get_domain_name(material_dict::MaterialDictType)::String
    mat_domains_master = get_domains_registry()
    if !haskey(material_dict, "mat_domain")
        return "GeneralDomain"
    end
    mat_domain = material_dict["mat_domain"]
    if !(mat_domain in mat_domains_master)
        throw(ArgumentError("Found invalid material domain type: $mat_domain"))
    end
    return mat_domain
end

function get_material_domains()::Vector{String}
    return get_domains_registry()
end

function get_type_name(material_dict::MaterialDictType)::String
    mat_types_master = get_types_registry()
    if !haskey(material_dict, "mat_type")
        return "General"
    end
    mat_type = material_dict["mat_type"]
    if !(mat_type in mat_types_master)
        throw(ArgumentError("Found invalid material type: $mat_type"))
    end
    return mat_type
end

function get_material_types()::Vector{String}
    return get_types_registry()
end

function check_and_clean_materials_input_dict(
    materials_dict::MaterialsDictType,
    material_option::Union{Nothing, Any} = nothing
)::MaterialsDictType
    checked_matids = String[]
    materials_dict_updated = MaterialsDictType()
    domains = String[]
    
    for (matid, material_dict) in materials_dict
        if !isa(material_dict, Dict)
            current_type = typeof(material_dict)
            throw(ArgumentError(
                "Material for ID $matid must be a dictionary. " *
                "The current type is $current_type. Adjust inputs."
            ))
        end
        
        matid = check_matid(matid)
        if matid in checked_matids
            throw(ArgumentError("Found duplicate material ID $matid. Adjust inputs."))
        end
        push!(checked_matids, matid)
        
        materials_dict_updated[matid] = check_and_clean_material(matid, material_dict)
        push!(domains, material_dict["mat_domain"])
    end
    add_layered_mantle_lithosphere_domains!(domains)
    check_required_domains(domains, material_option)
    return materials_dict_updated
end

function add_layered_mantle_lithosphere_domains!(domains::Vector{String})
    mdr = MaterialDomainsRegistry()
    if mdr.mantle_lithosphere in domains
        push!(domains, mdr.upper_mantle_lithosphere)
        push!(domains, mdr.middle_mantle_lithosphere)
        push!(domains, mdr.lower_mantle_lithosphere)
        print_info("Upper, middle, and lower mantle lithosphere domains have been added to the domains list.", level=1)
    end
    return nothing
end

function check_matid(matid::Union{String, Int16})::String
    if !(isa(matid, Int16) || isa(matid, String))
        throw(ArgumentError(
            "Material ID $matid is not a string or integer. Adjust inputs."
        ))
    end
    
    try
        parse(Int, string(matid))
    catch e
        if isa(e, ArgumentError)
            throw(ArgumentError(
                "Invalid material ID $matid. " *
                "Material ID could not be converted to an integer."
            ))
        end
        rethrow(e)
    end
    
    return string(matid)
end

function check_and_clean_material(
    matid::Union{String, Int16},
    material_dict::MaterialDictType
)::Dict{String, Union{Float64, Int64, String}}
    check_material_id(matid, material_dict)
    material_updated = check_parameter_names_and_types(material_dict)
    check_domain_and_material_type(material_dict)
    check_library_name(material_dict)
    return material_updated
end

function check_material_id(
    matid::Union{String, Int16},
    material_dict::MaterialDictType
)::Nothing
    if !haskey(material_dict, "matid")
        material_dict["matid"] = matid
    end
    return nothing
end

function check_parameter_names_and_types(
    material_dict::MaterialDictType
)::Dict{String, Union{Float64, Int64, String}}
    parameters_master = Dict{String, Type}(
        "matid" => String,
        "mat_name" => String,
        "mat_type" => String,
        "mat_domain" => String,
        "red_fraction" => Float64,
        "green_fraction" => Float64,
        "blue_fraction" => Float64
    )
    
    color_names = ["red_fraction", "green_fraction", "blue_fraction"]
    material_updated = Dict{String, Union{Float64, Int64, String}}()
    
    for (param_name, value) in material_dict
        if !haskey(parameters_master, param_name)
            throw(ArgumentError("$param_name is not a valid parameter name."))
        end
        material_updated[param_name] = convert(parameters_master[param_name], value)
    end
    
    for (param_name, _) in parameters_master
        if !haskey(material_updated, param_name)
            if param_name in color_names
                material_updated[param_name] = 0.5
            else
                throw(ArgumentError("material is missing required parameter $param_name"))
            end
        end
    end
    
    return material_updated
end

function check_domain_and_material_type(material_dict::MaterialDictType)::Nothing
    domains_master = get_domains_registry()
    types_master = get_types_registry()
    
    mat_domain = material_dict["mat_domain"]
    mat_type = material_dict["mat_type"]
    
    if !(mat_domain in domains_master)
        throw(ArgumentError("Invalid domain name: $mat_domain."))
    end
    if !(mat_type in types_master)
        throw(ArgumentError("Invalid type name: $mat_type."))
    end
    
    return nothing
end

function check_required_domains(
    domains::Vector{String},
    material_option::Union{Nothing, Any} = nothing
)::Nothing
    if material_option !== nothing && hasproperty(material_option, :required_domains)
        required_domains = material_option.required_domains
        if required_domains !== nothing
            missing_domains = String[]
            for required_domain in required_domains
                if !(required_domain in domains)
                    push!(missing_domains, required_domain)
                end
            end
            
            if !isempty(missing_domains)
                print_info("The following required domains were missing from materials_input_dict:", level=1)
                for missing_domain in missing_domains
                    print_info("$missing_domain", level=2)
                end
                print_info("", level=1)
                throw(ArgumentError("Missing required material domains. Adjust input."))
            end
        end
    end
    return nothing
end

function check_library_name(material_dict::MaterialDictType)::Nothing
    library_material_names_master = String[]
    mat_name = material_dict["mat_name"]
    if !(mat_name in library_material_names_master)
        # Note: In Python code this was just a pass statement
        # We could add validation here if needed
    end
    return nothing
end

end # module