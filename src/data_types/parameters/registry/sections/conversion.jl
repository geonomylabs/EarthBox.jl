function get_conversion_parameters()::NamedTuple
    return @params (
        KtoC = ParameterFloat(
            273.0, "None", "Conversion factor for Kelvins to Celsius: T_C = T_K - KtoC"),
        sec_per_yr = ParameterFloat(
            365.25*24.0*3600.0, "s/yr", "Seconds per year"),
        sec_per_Myr = ParameterFloat(
            365.25*24.0*3600.0*1e+6, "s/Myr", "Seconds per million years"),
        cm_yr2m_s = ParameterFloat(
            1.0/(100.0*365.25*24.0*3600.0), "m/s/cm/yr",
            "Conversion factor for cm/yr to m/s"
            ),
    )
end
