module MeltingTypes

function check_required_types_and_domains_for_melt_model(model::ModelData)::Nothing
    matid_types = model.materials.dicts.matid_types
    type_dict = Dict(
        "StickyAir" => matid_types["StickyAir"][1],
        "StickyWater" => matid_types["StickyWater"][1]
    )
    check_domain_dict(type_dict, "type", "contractional")

    domains = model.materials.dicts.matid_domains
    domain_dict = Dict(
        "Asthenosphere" => domains["Asthenosphere"],
        "SedimentaryBasin" => domains["SedimentaryBasin"],
        "BasalticOceanicCrust" => domains["BasalticOceanicCrust"],
        "GabbroicOceanicCrust" => domains["GabbroicOceanicCrust"],
        "MantleLithosphere" => domains["MantleLithosphere"],
        "WeakMantleLithosphere" => domains["WeakMantleLithosphere"]
    )
    check_domain_dict(domain_dict, "domain", "contractional")
    return nothing
end


end