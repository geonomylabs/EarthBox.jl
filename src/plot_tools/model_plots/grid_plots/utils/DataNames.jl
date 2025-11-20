module DataNames

Base.@kwdef struct ScalarH5Datanames
    # The following names are datanames from scalar field h5 files
    TempC::String = "TempC"
    log_eta_Pas::String = "log_eta_Pas"
    Eii_log::String = "Eii_log"
    pressure_GPa::String = "pressure_GPa"
    Sxx_MPa::String = "Sxx_MPa"
    Sxy_MPa::String = "Sxy_MPa"
    plastics::String = "plastics"
    plasticn::String = "plasticn"
    therm_cond::String = "therm_cond"
    rho_kg_m3::String = "rho_kg_m3"
    # The following names are datanames from velocity h5 files
    # Since they are individual scalar components, they are included with
    # the scalar datanames.
    velx_cmyr::String = "velx_cmyr"
    vely_cmyr::String = "vely_cmyr"
    velmag_cmyr::String = "velmag_cmyr"
end

function get_velocity_dataname_list(datanames::ScalarH5Datanames)
    return [datanames.velx_cmyr, datanames.vely_cmyr, datanames.velmag_cmyr]
end

end # module 