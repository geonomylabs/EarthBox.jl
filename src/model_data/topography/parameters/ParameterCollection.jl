module ParameterCollection

include("groups/__init__.jl")

import EarthBox.EarthBoxDtypes: AbstractParameterCollection
import .ModelOptionsGroup: ModelOptions
import .TopoGridGroup: TopoGrid
import .DepoAndErosionRatesGroup: DepoAndErosionRates
import .DownhillDiffusionGroup: DownhillDiffusion
import .SeaLevelGroup: SeaLevel
import .ReferenceLithosphereGroup: ReferenceLithosphere

"""
    Parameters <: AbstractParameterCollection

Parameter collection for topography evolution model.

# Fields
- `model_options::`[`ModelOptions`](@ref): Topography model on/off switches and options
- `topo_grid::`[`TopoGrid`](@ref): Topography grid parameters
- `depo_and_erosion_rates::`[`DepoAndErosionRates`](@ref): Deposition and erosion rates
- `downhill_diffusion::`[`DownhillDiffusion`](@ref): Downhill diffusion parameters
- `sealevel::`[`SeaLevel`](@ref): Sea level parameters
- `reference_lithosphere::`[`ReferenceLithosphere`](@ref): Reference lithosphere parameters

# Constructor
    Parameters()

Create a new Parameters collection with default topography model values.

"""
mutable struct Parameters <: AbstractParameterCollection
    model_options::ModelOptions
    topo_grid::TopoGrid
    depo_and_erosion_rates::DepoAndErosionRates
    downhill_diffusion::DownhillDiffusion
    sealevel::SeaLevel
    reference_lithosphere::ReferenceLithosphere
end

function Parameters()::Parameters
    return Parameters(
        ModelOptions(),
        TopoGrid(),
        DepoAndErosionRates(),
        DownhillDiffusion(),
        SeaLevel(),
        ReferenceLithosphere()
    )
end

end # module 