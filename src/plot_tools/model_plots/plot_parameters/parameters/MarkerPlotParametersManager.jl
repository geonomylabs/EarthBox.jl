module MarkerPlotParametersManager

import ...PlotDtypes: PlotDictType
import ...PlotDtypes: AbstractPlotParameterGroup

""" Marker plot parameters class for setting marker plot parameters.

Attributes are initialized with default values that are read from a
template stored in plot_dict. The template is defined in the following
yamal file:

    plot_tools/model_plots/utils/plot_dict_template.yml

"""
Base.@kwdef mutable struct MarkerPlotParameters <: AbstractPlotParameterGroup
    # Miscellaneous parameters
    colorbar_shift_factor::Float64 = 1.0
    marker_size::Float64 = 0.0
    plot_contour_labels::Int = 0
    contour_line_width::Float64 = 0.0
    contour_line_color::String = "black"
    decimation_factor_scatter_overlay::Int = 0
    nx_contour_grid::Int = 0
    ny_contour_grid::Int = 0
    decimation_factor::Int = 1
    # Basic grid mesh
    plot_mesh::Int = 0
    mesh_line_width::Float64 = 0.0
    # Topography marker chain
    plot_topography::Int = 0
    topo_line_width::Float64 = 0.0
    topo_line_color::String = "black"
    # Base level plot parameters
    plot_base_level::Int = 0
    base_level_line_width::Float64 = 0.0
    base_level_line_color::String = "blue"
    # Marker plastic strain
    plot_plastic_strain::Int = 0
    plot_plastic_strain_contours::Int = 0
    strain_min::Float64 = 0.0
    strain_max::Float64 = 0.0
    strain_contour_interval::Float64 = 0.0
    strain_cmap::String = "None"
    # Marker plastic strain sate
    plot_plastic_strain_rate::Int = 0
    plot_plastic_strain_rate_contours::Int = 0
    strain_rate_min::Float64 = 0.0
    strain_rate_max::Float64 = 0.0
    strain_rate_contour_interval::Float64 = 0.0
    strain_rate_cmap::String = "None"
    # Marker sediment age
    plot_sediment_age::Int = 0
    plot_sediment_age_contours::Int = 0
    age_min::Float64 = 0.0
    age_max::Float64 = 0.0
    age_contour_interval::Float64 = 0.0
    age_cmap::String = "None"
    plot_volcanics_age::Int = 0
    plot_volcanics_age_contours::Int = 0
    age_min_volcanics::Float64 = 0.0
    age_max_volcanics::Float64 = 0.0
    age_contour_interval_volcanics::Float64 = 0.0
    use_alternating_colormap_volcanics::Bool = true
    age_cmap_volcanics::String = "None"
    plot_intrusive_age::Int = 0
    plot_intrusive_age_contours::Int = 0
    age_min_intrusive::Float64 = 0.0
    age_max_intrusive::Float64 = 0.0
    age_contour_interval_intrusive::Float64 = 0.0
    use_alternating_colormap_intrusive::Bool = true
    age_cmap_intrusive::String = "None"
    # Temperature contours on basic grid
    plot_temperature_contours::Int = 0
    temperature_min::Float64 = 0.0
    temperature_max::Float64 = 0.0
    temperature_contour_interval::Float64 = 0.0
    temperature_contour_color::String = "black"
    temperature_number_format::String = "%6.2f"
    temperature_label_rightside_up::Bool = true
    # Marker melt fraction
    plot_meltfrac::Int = 0
    plot_meltfrac_contours::Int = 0
    # Marker extracted melt fraction
    plot_extracted_meltfrac::Int = 0
    plot_extracted_meltfrac_contours::Int = 0
    # Marker extractable melt fraction
    plot_extractable_meltfrac::Int = 0
    plot_extractable_meltfrac_contours::Int = 0
    # Parameters that apply to all meltfrac related marker plots
    melt_fraction_min::Float64 = 0.0
    melt_fraction_max::Float64 = 0.0
    melt_fraction_contour_interval::Float64 = 0.0
    meltfrac_contour_color::String = "red"
    meltfrac_cmap::String = "None"
    meltfrac_number_format::String = "%6.2f"
    meltfrac_label_rightside_up::Bool = true
    # Marker serpentinization
    plot_serpentinization::Int = 0
    plot_serpentinization_contours::Int = 0
    serpentinization_min::Float64 = 0.0
    serpentinization_max::Float64 = 0.0
    serpentinization_contour_interval::Float64 = 0.0
    serpentinization_cmap::String = "None"
    # Marker density
    plot_density::Int = 0
    plot_density_contours::Int = 0
    # Marker mantle density
    plot_mantle_density::Int = 0
    plot_mantle_density_contours::Int = 0
    # Marker density parameters used for all marker density plots
    density_min::Float64 = 0.0
    density_max::Float64 = 0.0
    density_contour_interval::Float64 = 0.0
    density_cmap::String = "None"
    # Heat flow subplot parameters
    heatflow_min::Float64 = 0.0
    heatflow_max::Float64 = 200.0
    heatflow_spacing::Float64 = 10.0
    # Gravity subplot parameters
    gravity_min::Float64 = -200.0
    gravity_max::Float64 = 200.0
    gravity_spacing::Float64 = 50.0
    # Subplot parameters
    height_ratios::Vector{Float64} = [0.25, 0.25, 0.75]
    subplot_spacing_fraction::Float64 = 0.4
end

function MarkerPlotParameters(plot_dict::PlotDictType)::MarkerPlotParameters
    plot_params = plot_dict["general_parameters"]
    return MarkerPlotParameters(
        colorbar_shift_factor = get(plot_params, "colorbar_shift_factor", 1.0),
        marker_size = get(plot_params, "marker_size", 0.0),
        plot_contour_labels = get(plot_params, "plot_contour_labels", 0),
        contour_line_width = get(plot_params, "contour_line_width", 0.0),
        contour_line_color = get(plot_params, "contour_line_color", "black"),
        decimation_factor_scatter_overlay = get(plot_params, "decimation_factor_scatter_overlay", 0),
        nx_contour_grid = get(plot_params, "nx_contour_grid", 0),
        ny_contour_grid = get(plot_params, "ny_contour_grid", 0),
        decimation_factor = get(plot_params, "decimation_factor", 1),
        plot_mesh = get(plot_params, "plot_mesh", 0),
        mesh_line_width = get(plot_params, "mesh_line_width", 0.0),
        plot_topography = get(plot_params, "plot_topography", 0),
        topo_line_width = get(plot_params, "topo_line_width", 0.0),
        topo_line_color = get(plot_params, "topo_line_color", "black"),
        plot_base_level = get(plot_params, "plot_base_level", 0),
        base_level_line_width = get(plot_params, "base_level_line_width", 0.0),
        base_level_line_color = get(plot_params, "base_level_line_color", "blue"),
        plot_plastic_strain = get(plot_params, "plot_plastic_strain", 0),
        plot_plastic_strain_contours = get(plot_params, "plot_plastic_strain_contours", 0),
        strain_min = get(plot_params, "strain_min", 0.0),
        strain_max = get(plot_params, "strain_max", 0.0),
        strain_contour_interval = get(plot_params, "strain_contour_interval", 0.0),
        strain_cmap = get(plot_params, "strain_cmap", "None"),
        plot_plastic_strain_rate = get(plot_params, "plot_plastic_strain_rate", 0),
        plot_plastic_strain_rate_contours = get(plot_params, "plot_plastic_strain_rate_contours", 0),
        strain_rate_min = get(plot_params, "strain_rate_min", 0.0),
        strain_rate_max = get(plot_params, "strain_rate_max", 0.0),
        strain_rate_contour_interval = get(plot_params, "strain_rate_contour_interval", 0.0),
        strain_rate_cmap = get(plot_params, "strain_rate_cmap", "None"),
        plot_sediment_age = get(plot_params, "plot_sediment_age", 0),
        plot_sediment_age_contours = get(plot_params, "plot_sediment_age_contours", 0),
        age_min = get(plot_params, "age_min", 0.0),
        age_max = get(plot_params, "age_max", 0.0),
        age_contour_interval = get(plot_params, "age_contour_interval", 0.0),
        age_cmap = get(plot_params, "age_cmap", "None"),
        plot_volcanics_age = get(plot_params, "plot_volcanics_age", 0),
        plot_volcanics_age_contours = get(plot_params, "plot_volcanics_age_contours", 0),
        age_min_volcanics = get(plot_params, "age_min_volcanics", 0.0),
        age_max_volcanics = get(plot_params, "age_max_volcanics", 0.0),
        age_contour_interval_volcanics = get(plot_params, "age_contour_interval_volcanics", 0.0),
        use_alternating_colormap_volcanics = get(plot_params, "use_alternating_colormap_volcanics", true),
        age_cmap_volcanics = get(plot_params, "age_cmap_volcanics", "None"),
        plot_intrusive_age = get(plot_params, "plot_intrusive_age", 0),
        plot_intrusive_age_contours = get(plot_params, "plot_intrusive_age_contours", 0),
        age_min_intrusive = get(plot_params, "age_min_intrusive", 0.0),
        age_max_intrusive = get(plot_params, "age_max_intrusive", 0.0),
        age_contour_interval_intrusive = get(plot_params, "age_contour_interval_intrusive", 0.0),
        use_alternating_colormap_intrusive = get(plot_params, "use_alternating_colormap_intrusive", true),
        age_cmap_intrusive = get(plot_params, "age_cmap_intrusive", "None"),
        plot_temperature_contours = get(plot_params, "plot_temperature_contours", 0),
        temperature_min = get(plot_params, "temperature_min", 0.0),
        temperature_max = get(plot_params, "temperature_max", 0.0),
        temperature_contour_interval = get(plot_params, "temperature_contour_interval", 0.0),
        temperature_contour_color = get(plot_params, "temperature_contour_color", "black"),
        temperature_number_format = get(plot_params, "temperature_number_format", "%6.2f"),
        temperature_label_rightside_up = get(plot_params, "temperature_label_rightside_up", true),
        plot_meltfrac = get(plot_params, "plot_meltfrac", 0),
        plot_meltfrac_contours = get(plot_params, "plot_meltfrac_contours", 0),
        plot_extracted_meltfrac = get(plot_params, "plot_extracted_meltfrac", 0),
        plot_extracted_meltfrac_contours = get(plot_params, "plot_extracted_meltfrac_contours", 0),
        plot_extractable_meltfrac = get(plot_params, "plot_extractable_meltfrac", 0),
        plot_extractable_meltfrac_contours = get(plot_params, "plot_extractable_meltfrac_contours", 0),
        melt_fraction_min = get(plot_params, "melt_fraction_min", 0.0),
        melt_fraction_max = get(plot_params, "melt_fraction_max", 0.0),
        melt_fraction_contour_interval = get(plot_params, "melt_fraction_contour_interval", 0.0),
        meltfrac_contour_color = get(plot_params, "meltfrac_contour_color", "red"),
        meltfrac_cmap = get(plot_params, "meltfrac_cmap", "None"),
        meltfrac_number_format = get(plot_params, "meltfrac_number_format", "%6.2f"),
        meltfrac_label_rightside_up = get(plot_params, "meltfrac_label_rightside_up", true),
        plot_serpentinization = get(plot_params, "plot_serpentinization", 0),
        plot_serpentinization_contours = get(plot_params, "plot_serpentinization_contours", 0),
        serpentinization_min = get(plot_params, "serpentinization_min", 0.0),
        serpentinization_max = get(plot_params, "serpentinization_max", 0.0),
        serpentinization_contour_interval = get(plot_params, "serpentinization_contour_interval", 0.0),
        serpentinization_cmap = get(plot_params, "serpentinization_cmap", "None"),
        plot_density = get(plot_params, "plot_density", 0),
        plot_density_contours = get(plot_params, "plot_density_contours", 0),
        density_min = get(plot_params, "density_min", 0.0),
        density_max = get(plot_params, "density_max", 0.0),
        density_contour_interval = get(plot_params, "density_contour_interval", 0.0),
        density_cmap = get(plot_params, "density_cmap", "None"),
        heatflow_min = get(plot_params, "heatflow_min", 0.0),
        heatflow_max = get(plot_params, "heatflow_max", 200.0),
        heatflow_spacing = get(plot_params, "heatflow_spacing", 10.0),
        gravity_min = get(plot_params, "gravity_min", -200.0),
        gravity_max = get(plot_params, "gravity_max", 200.0),
        gravity_spacing = get(plot_params, "gravity_spacing", 50.0),
        height_ratios = get(plot_params, "height_ratios", [0.25, 0.25, 0.75]),
        subplot_spacing_fraction = get(plot_params, "subplot_spacing_fraction", 0.4),
    )
end

end # module