module ConversionFuncs

"""
    meters_squared_per_year_to_meters_squared_per_second(m2_yr::Float64)::Float64

Convert m²/yr to m²/s.
"""
function meters_squared_per_year_to_meters_squared_per_second(m2_yr::Float64)::Float64
    m2_s = m2_yr/(365.25*24.0*3600.0)
    return m2_s
end

function m2_yr_to_m2_s(m2_yr::Float64)::Float64
    return meters_squared_per_year_to_meters_squared_per_second(m2_yr)
end

"""
    meters_squared_per_second_to_meters_squared_per_year(m2_s::Float64)::Float64

Convert m²/s to m²/yr.
"""
function meters_squared_per_second_to_meters_squared_per_year(m2_s::Float64)::Float64
    m2_yr = m2_s*(365.25*24.0*3600.0)
    return m2_yr
end

function m2_s_to_m2_yr(m2_s::Float64)::Float64
    return meters_squared_per_second_to_meters_squared_per_year(m2_s)
end

"""
    millions_of_years_to_seconds(time_myr::Float64)::Float64

Convert millions of years to seconds.
"""
function millions_of_years_to_seconds(time_myr::Float64)::Float64
    time_years = time_myr*1e6
    time_seconds = years_to_seconds(time_years)
    return time_seconds
end

"""
    years_to_seconds(time_years::Float64)::Float64

Convert years to seconds.
"""
@inline function years_to_seconds(time_years::Float64)::Float64
    seconds_per_year = 365.25*24.0*3600.0
    time_seconds = time_years*seconds_per_year
    return time_seconds
end

function yr_to_s(time_years::Float64)::Float64
    return years_to_seconds(time_years)
end

"""
    seconds_to_years(time_seconds::Float64)::Float64

Convert seconds to years.
"""
@inline function seconds_to_years(time_seconds::Float64)::Float64
    seconds_per_year = 365.25*24.0*3600.0
    time_years = time_seconds/seconds_per_year
    return time_years
end

"""
    celsius_to_kelvin(temperature_celsius::Float64)::Float64

Convert Celsius to Kelvins.
"""
@inline function celsius_to_kelvin(temperature_celsius::Float64)::Float64
    temperature_kelvins = temperature_celsius + 273.0
    return temperature_kelvins
end

"""
    kelvin_to_celsius(temperature_kelvins::Float64)::Float64

Convert Kelvins to Celsius.
"""
@inline function kelvin_to_celsius(temperature_kelvins::Float64)::Float64
    temperature_celsius = temperature_kelvins - 273.0
    return temperature_celsius
end

"""
    get_factor_cm_yr_to_m_s()::Float64

Get factor to convert cm/yr to m/s.
"""
function get_factor_cm_yr_to_m_s()::Float64
    return 1.0/100.0/(365.25*24.0*3600.0)
end

"""
    cm_yr_to_m_s(velext_cm_yr::Float64)::Float64

Convert cm/yr to m/s.
"""
function cm_yr_to_m_s(velext_cm_yr::Float64)::Float64
    return cm_per_yr_to_meters_per_seconds(velext_cm_yr)
end

"""
    cm_per_yr_to_meters_per_seconds(velext_cm_yr::Float64)::Float64

Convert cm/yr to m/s.
"""
function cm_per_yr_to_meters_per_seconds(velext_cm_yr::Float64)::Float64
    cm_yr2m_s = get_factor_cm_yr_to_m_s()
    velext_m_s = velext_cm_yr*cm_yr2m_s
    return velext_m_s
end

"""
    m_s_to_cm_yr(velext_m_s::Float64)::Float64

Convert m/s to cm/yr.
"""
function m_s_to_cm_yr(velext_m_s::Float64)::Float64
    return meters_per_seconds_to_cm_per_yr(velext_m_s)
end

"""
    meters_per_seconds_to_cm_per_yr(velext_m_s::Float64)::Float64

Convert m/s to cm/yr.
"""
function meters_per_seconds_to_cm_per_yr(velext_m_s::Float64)::Float64
    cm_yr2m_s = get_factor_cm_yr_to_m_s()
    velext_cm_yr = velext_m_s/cm_yr2m_s
    return velext_cm_yr
end

"""
    meters_per_seconds_to_mm_per_yr(velext_m_s::Float64)::Float64

Convert m/s to mm/yr.
"""
function meters_per_seconds_to_mm_per_yr(velext_m_s::Float64)::Float64
    cm_yr2m_s = get_factor_cm_yr_to_m_s()
    velext_mm_yr = velext_m_s/cm_yr2m_s*10.0
    return velext_mm_yr
end

"""
    meters_per_year_to_meters_per_seconds(rate_m_yr::Float64)::Float64

Convert m/yr to m/s.
"""
function meters_per_year_to_meters_per_seconds(rate_m_yr::Float64)::Float64
    return rate_m_yr/(365.25*24.0*3600.0)
end

function m_yr_to_m_s(rate_m_yr::Float64)::Float64
    return meters_per_year_to_meters_per_seconds(rate_m_yr)
end

"""
    mm_per_yr_to_meters_per_seconds(rate_mm_yr::Float64)::Float64

Convert mm/yr to m/s.
"""
function mm_per_yr_to_meters_per_seconds(rate_mm_yr::Float64)::Float64
    cm_yr2m_s = get_factor_cm_yr_to_m_s()
    rate_m_s = rate_mm_yr*cm_yr2m_s/10.0
    return rate_m_s
end

function mm_yr_to_m_s(rate_mm_yr::Float64)::Float64
    return mm_per_yr_to_meters_per_seconds(rate_mm_yr)
end

"""
    meters_per_second_squared_to_gal(acceleration_m_s2::Float64)::Float64

Convert m/s² to gal (cm/s²).
"""
function meters_per_second_squared_to_gal(acceleration_m_s2::Float64)::Float64
    acceleration_gal = acceleration_m_s2*1e2
    return acceleration_gal
end

"""
    gal_to_milligal(acceleration_gal::Float64)::Float64

Convert gal to mGal.
"""
function gal_to_milligal(acceleration_gal::Float64)::Float64
    acceleration_mgal = acceleration_gal*1e+3
    return acceleration_mgal
end

"""
    meters_per_second_squared_to_milligal(acceleration_m_s2::Float64)::Float64

Convert m/s² to milligal (mGal).
"""
function meters_per_second_squared_to_milligal(acceleration_m_s2::Float64)::Float64
    acceleration_gal = meters_per_second_squared_to_gal(acceleration_m_s2)
    acceleration_mgal = gal_to_milligal(acceleration_gal)
    return acceleration_mgal
end

export meters_squared_per_year_to_meters_squared_per_second,
       meters_squared_per_second_to_meters_squared_per_year,
       millions_of_years_to_seconds,
       years_to_seconds,
       seconds_to_years,
       celsius_to_kelvin,
       kelvin_to_celsius,
       cm_per_yr_to_meters_per_seconds,
       get_factor_cm_yr_to_m_s,
       meters_per_seconds_to_cm_per_yr,
       meters_per_seconds_to_mm_per_yr,
       meters_per_year_to_meters_per_seconds,
       mm_per_yr_to_meters_per_seconds,
       meters_per_second_squared_to_gal,
       gal_to_milligal,
       meters_per_second_squared_to_milligal

end # module 