module ParameterCollection

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import EarthBox.Parameters: ParameterFloat

const ROOT_NAME = "model.conversion.parameters"

const PDATA = get_eb_parameters()

"""
    Parameters <: AbstractParameterCollection

Parameter collection for conversion parameters.

# Fields
- `KtoC::`[`ParameterFloat`](@ref): $(PDATA.KtoC.description)
- `sec_per_yr::`[`ParameterFloat`](@ref): $(PDATA.sec_per_yr.description)
- `sec_per_Myr::`[`ParameterFloat`](@ref): $(PDATA.sec_per_Myr.description)
- `cm_yr2m_s::`[`ParameterFloat`](@ref): $(PDATA.cm_yr2m_s.description)

# Constructor
    Parameters()

# Returns
- `Parameters`: New Parameters collection with initialized values

# Nested Dot Access
- `KtoC = $(ROOT_NAME).KtoC.value`
- `sec_per_yr = $(ROOT_NAME).sec_per_yr.value`
- `sec_per_Myr = $(ROOT_NAME).sec_per_Myr.value`
- `cm_yr2m_s = $(ROOT_NAME).cm_yr2m_s.value`

"""
mutable struct Parameters <: AbstractParameterCollection
    KtoC::ParameterFloat
    sec_per_yr::ParameterFloat
    sec_per_Myr::ParameterFloat
    cm_yr2m_s::ParameterFloat
end

function Parameters()::Parameters
    pdata = get_eb_parameters()
    return Parameters(
        pdata.KtoC,
        pdata.sec_per_yr,
        pdata.sec_per_Myr,
        pdata.cm_yr2m_s
    )
end

end # module
