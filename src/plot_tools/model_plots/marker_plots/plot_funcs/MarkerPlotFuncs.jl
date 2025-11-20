module MarkerPlotFuncs

include("plots/CompositionPlot.jl")
include("plots/Serpentinization.jl")
include("plots/MarkerAge.jl")
include("plots/Strain.jl")
include("plots/PlasticStrainRate.jl")
include("plots/ExtractableMeltfrac.jl")
include("plots/ExtractedMeltfrac.jl")
include("plots/Mesh.jl")
include("plots/TopographyPlot.jl")
include("plots/TemperatureContours.jl")
include("plots/PlasticFailure.jl")
include("plots/Meltfrac.jl")
include("plots/Density.jl")

import EarthBox.Markers.MarkerMaterials.MaterialsContainer: Materials
import .CompositionPlot: plot_composition
import .TemperatureContours
import .TopographyPlot
import ..PlotMarkerArraysManager: PlotMarkerArrays
import ..PlotTopoArraysManager: PlotTopoArrays
import ..MarkerColormapManager: MarkerColorMap
import ...PlotParametersManager: PlotParameters
import ...PlotParametersManager: reset_plot_counter!
import ...PlotParametersManager.PlotContoursManager: reset_contour_description!
import ...GridPlotsManager.PlotScalarArraysManager: PlotScalarArrays
import ...PlotDtypes: AxesType

function plot_marker_scalars(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    materials::Materials,
    colormap::MarkerColorMap,
    axes::AxesType
)::Nothing
    reset_plot_counter!(parameters)
    reset_contour_description!(parameters.contours)
    plot_composition(parameters, marker_arrays, materials, colormap, axes)
    MarkerAge.plot_filtered_sediment_age(parameters, marker_arrays, materials, axes)
    MarkerAge.plot_filtered_volcanics_age(parameters, marker_arrays, materials, axes)
    MarkerAge.plot_filtered_intrusive_age(parameters, marker_arrays, materials, axes)
    Strain.plot_filtered_strain(parameters, marker_arrays, materials, axes)
    PlasticStrainRate.plot_filtered_plastic_strain_rate(parameters, marker_arrays, axes)
    Serpentinization.plot_filtered_serpentinization(parameters, marker_arrays, axes)
    Density.plot_filtered_mantle_density(parameters, marker_arrays, materials, axes)
    Meltfrac.plot_filtered_meltfrac(parameters, marker_arrays, materials, axes)
    ExtractedMeltfrac.plot_filtered_extracted_meltfrac(parameters, marker_arrays, materials, axes)
    ExtractableMeltfrac.plot_filtered_extractable_meltfrac(parameters, marker_arrays, materials, axes)
    return nothing
end

function count_plots!(parameters::PlotParameters)::Nothing
    total_plots =  (
        1 # composition is always plotted
        + parameters.marker_plot_params.plot_sediment_age # sediment age
        + parameters.marker_plot_params.plot_volcanics_age # volcanics age
        + parameters.marker_plot_params.plot_intrusive_age # intrusive age
        + parameters.marker_plot_params.plot_plastic_strain # plastic strain
        + parameters.marker_plot_params.plot_plastic_strain_rate # plastic strain rate
        + parameters.marker_plot_params.plot_serpentinization # serpentinization
        + parameters.marker_plot_params.plot_mantle_density # mantle density
        + parameters.marker_plot_params.plot_meltfrac # melt fraction
        + parameters.marker_plot_params.plot_extracted_meltfrac # extracted melt fraction
        + parameters.marker_plot_params.plot_extractable_meltfrac # extractable melt fraction
    )
    parameters.total_plots = total_plots
    return nothing
end

function plot_density(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    materials::Materials,
    axes::AxesType
)::Nothing
    reset_plot_counter!(parameters)
    reset_contour_description!(parameters.contours)
    Density.plot_filtered_density(parameters, marker_arrays, axes)
    Meltfrac.plot_filtered_meltfrac(parameters, marker_arrays, materials, axes)
    return nothing
end

function plot_mesh(
    parameters::PlotParameters,
    scalar_arrays::PlotScalarArrays,
    axes::AxesType,
)::Nothing
    Mesh.plot_mesh(parameters, scalar_arrays, axes)
    return nothing
end

function plot_plastic_failure(
    parameters::PlotParameters,
    marker_arrays::PlotMarkerArrays,
    colormap::Any,
    axes::AxesType
)::Nothing
    PlasticFailure.plot_plastic(parameters, marker_arrays, colormap, axes)
    return nothing
end

function plot_topo(
    parameters::PlotParameters,
    topo_arrays::PlotTopoArrays,
    axes::AxesType
)::Nothing
    TopographyPlot.plot_topo(parameters, topo_arrays, axes)
    return nothing
end

function plot_base_level(
    parameters::PlotParameters,
    topo_arrays::PlotTopoArrays,
    axes::AxesType
)::Nothing
    TopographyPlot.plot_base_level(parameters, topo_arrays, axes)
    return nothing
end
 
function plot_temperature_contours(
    parameters::PlotParameters,
    scalar_arrays::PlotScalarArrays,
    axes::AxesType
)::Nothing
    TemperatureContours.plot_temperature_contours!(parameters, scalar_arrays, axes)
    return nothing
end

end # module
