module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import ...Grids2dContainer: Grids
import .DistributionGroup: Distribution
import .AdvectionGroup: Advection
import .SubgridDiffusionGroup: SubgridDiffusion
import .RecyclingGroup: Recycling

"""
    Parameters <: AbstractParameterCollection

Collection of marker parameters.

# Fields
- `distribution::`[`Distribution`](@ref): Marker distribution parameters
- `advection::`[`Advection`](@ref): Marker advection parameters
- `subgrid_diffusion::`[`SubgridDiffusion`](@ref): Subgrid diffusion parameters
- `recycling::`[`Recycling`](@ref): Marker recycling parameters
"""
mutable struct Parameters <: AbstractParameterCollection
    distribution::Distribution
    advection::Advection
    subgrid_diffusion::SubgridDiffusion
    recycling::Recycling
end

function Parameters(
    marker_parameters::NamedTuple
)::Parameters
    return Parameters(
        Distribution(marker_parameters),
        Advection(),
        SubgridDiffusion(),
        Recycling()
    )
end

end # module 