module ParameterCollection

import EarthBox.ParameterRegistry: get_eb_parameters
import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import EarthBox.Parameters: ParameterInt

const ROOT_NAME = "model.benchmarks.parameters"

const PDATA = get_eb_parameters()

"""
    Parameters <: AbstractParameterCollection

Parameter collection for benchmark parameters.

# Fields
- `iuse_ramberg_post_processing::`[`ParameterInt`](@ref): $(PDATA.iuse_ramberg_post_processing.description)
- `iuse_viscous_block_processing::`[`ParameterInt`](@ref): $(PDATA.iuse_viscous_block_processing.description)
- `iuse_conbox_post_processing::`[`ParameterInt`](@ref): $(PDATA.iuse_conbox_post_processing.description)

# Nested Dot Access
- `iuse_ramberg_post_processing = $(ROOT_NAME).iuse_ramberg_post_processing.value`
- `iuse_viscous_block_processing = $(ROOT_NAME).iuse_viscous_block_processing.value`
- `iuse_conbox_post_processing = $(ROOT_NAME).iuse_conbox_post_processing.value`

# Constructor
    Parameters()

# Returns
- `Parameters`: New Parameters collection with initialized values
"""
mutable struct Parameters <: AbstractParameterCollection
    iuse_ramberg_post_processing::ParameterInt
    iuse_viscous_block_processing::ParameterInt
    iuse_conbox_post_processing::ParameterInt
end

function Parameters()::Parameters
    pdata = get_eb_parameters()
    return Parameters(
        pdata.iuse_ramberg_post_processing,
        pdata.iuse_viscous_block_processing,
        pdata.iuse_conbox_post_processing
    )
end

end # module
