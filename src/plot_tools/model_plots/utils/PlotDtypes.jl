module PlotDtypes

import CairoMakie

abstract type AbstractPlotParameterGroup end

const AxesType = CairoMakie.Axis

const PlotParametersType = Dict{String, Union{Float64, Int64, String, Bool, Tuple{Float64, Vararg{Float64}}}}
const PlotDictType = Dict{String, PlotParametersType}
const ParameterListType = Vector{Union{Float64, Int64, String, Bool, Tuple{Float64, Vararg{Float64}}, String}}
const PlotTemplateParametersType = Dict{String, ParameterListType}
const PlotDictTemplateType = Dict{String, PlotTemplateParametersType}


end # module 