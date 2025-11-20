module LiquidusModels

import EarthBox.OptionTools: OptionState
import EarthBox.OptionTools: make_option_names

const option_ids = Dict{Symbol, Int}(
    :NoMelting => -1,
    :SedimentGerya2010 => 0,
    :UprContCrustGerya2010 => 1,
    :LwrContCrustGerya2010 => 2,
    :DryMantleGerya2010 => 3,
    :PeridotiteKatz2003 => 4,
    :GabbroGerya2010 => 5,
    :GabbroGlacier => 6
)

const option_names = make_option_names(option_ids)

function get_options()::Dict{Int, OptionState}
    Dict{Int, OptionState}(
        option_ids[option_names.NoMelting] => 
            OptionState(
                option_name=string(option_names.NoMelting),
                description="Negative one is used to signal a material without melting. " *
                           "This is the default value for itype_liquidus."
            ),
        option_ids[option_names.SedimentGerya2010] => 
            OptionState(
                option_name=string(option_names.SedimentGerya2010),
                description="Sediment liquidus model from [gerya2010](@citet) Table 17.2."
            ),
        option_ids[option_names.UprContCrustGerya2010] => 
            OptionState(
                option_name=string(option_names.UprContCrustGerya2010),
                description="Upper continental crust liquidus model from [gerya2010](@citet) Table 17.2."
            ),
        option_ids[option_names.LwrContCrustGerya2010] => 
            OptionState(
                option_name=string(option_names.LwrContCrustGerya2010),
                description="Lower continental crust liquidus model from [gerya2010](@citet) Table 17.2."
            ),
        option_ids[option_names.DryMantleGerya2010] => 
            OptionState(
                option_name=string(option_names.DryMantleGerya2010),
                description="Dry mantle liquidus model from [gerya2010](@citet) Table 17.2."
            ),
        option_ids[option_names.PeridotiteKatz2003] => 
            OptionState(
                option_name=string(option_names.PeridotiteKatz2003),
                description="Peridotite liquidus model from [katz03](@citet)."
            ),
        option_ids[option_names.GabbroGerya2010] => 
            OptionState(
                option_name=string(option_names.GabbroGerya2010),
                description="Gabbro liquidus model from [gerya2010](@citet) Table 17.2."
            ),
        option_ids[option_names.GabbroGlacier] => 
            OptionState(
                option_name=string(option_names.GabbroGlacier),
                description="Gabbro liquidus model that approximates the gabbro glacier " *
                           "model of [maclennan04](@citet)."
            )
    )
end

function format_liquidus_options()
    liquidus_options = get_options()
    sorted_ids = sort(collect(keys(liquidus_options)))
    
    # Create header
    header = "| `stype_liquidus` | `itype_liquidus` | Description |\n"
    separator = "|:-----------------|:----------------:|:------------|\n"
    
    # Create rows
    rows = join(
        ["| `\"$(liquidus_options[id].option_name)\"` | $(id) | " *
         "$(liquidus_options[id].description) |" 
         for id in sorted_ids], 
        "\n"
    )
    
    return header * separator * rows
end

end # module 