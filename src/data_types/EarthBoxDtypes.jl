module EarthBoxDtypes

import DataStructures: OrderedDict

abstract type AbstractXDMFTimeStep end
abstract type AbstractOptionType end
abstract type AbstractEarthBoxArray1D end
abstract type AbstractEarthBoxArray2D end
abstract type AbstractEarthBoxArray3D end
abstract type AbstractArrayGroup end
abstract type AbstractArrayCollection end
abstract type AbstractParameter end
abstract type AbstractParameterGroup end
abstract type AbstractParameterCollection end
abstract type CollectionContainer end

const InputParamDictType = Dict{String, Tuple{Function, Union{Float64, Int64, String}, String}}
const InputMaterialParamDictType = Dict{String, Tuple{Function, Any, String}}
const ParametersDictType = Dict{String, Vector{Any}}

"""
    MaterialDictType

Used to store a collection of material properties for a single material where the key is the 
property name and the value is the property value.
"""
const MaterialDictType = Dict{String, Union{Float64, Int16, Int64, String}}

"""
    MaterialsDictType

Used to store a collection of materials where the key is the material name or integer ID and the
value is a dictionary of material properties where the key is the property name and the value is the 
property value.
"""
const MaterialsDictType = Dict{Union{String, Int16}, MaterialDictType}

const RawMaterialParameters = Dict{String, Vector{Union{Float64, Int64, String}}}
const RawMaterialsInputType = Dict{String, Union{String, RawMaterialParameters}}
const MasterMaterialParametersType = Dict{String, Vector{Union{Function, Float64, Int64, String}}}

# Type for material model dictionary composed of multiple materials
const MaterialModelDictType = Dict{String, MaterialDictType}

# Type for dictionary storing input for a single material
const MaterialInputDictType = Dict{String, Vector{Union{Float64, Int64, String}}}
# Type for full dictionary of input material information read
# directly from yaml file
const MaterialsInputDictType = Dict{Union{Int16, String}, Union{MaterialInputDictType, String, Int16}}
# Type for material library read from yaml file
const MaterialsLibraryDictType = Dict{Union{Int16, String}, Union{MaterialModelDictType, String, Int16}}
# Type for a section of model input file
const ModelInputSectionDictType = Dict{String, Vector{Union{Float64, Int64, String}}}
# Type for full model input dictionary read from yaml file
const ModelInputDictType = Dict{String, ModelInputSectionDictType}

# Type for test information
const TestInfoDictType = Dict{String, Dict{String, Union{Vector{Int64}, Int64, String}}}

# Types used for plotting
const PlotParametersType = Dict{String, Vector{Union{Float64, Int64, String}}}
const PlotDictType = Dict{String, PlotParametersType}

# Type for initialization parameters
const InitParamDictType = Dict{String, Vector{Union{Float64, Int64, String}}}

const InputDictType = Dict{String, Union{Float64, Int64, String}}

const ObjDictType = Dict{String, Union{AbstractParameter, AbstractEarthBoxArray1D, AbstractEarthBoxArray2D}}

# Picard loop information
const GlobalIterDictType = OrderedDict{Any, Any}

end # module 