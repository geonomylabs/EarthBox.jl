module PlotViewManager

import ...PlotDtypes: PlotDictType, AbstractPlotParameterGroup

Base.@kwdef mutable struct PlotView <: AbstractPlotParameterGroup
    use_close_up::Union{Int, Bool}
    xmin_normal::Float64
    xmax_normal::Float64
    ymin_normal::Float64
    ymax_normal::Float64
    xmin_close_up::Float64
    xmax_close_up::Float64
    ymin_close_up::Float64
    ymax_close_up::Float64
    xmin_active::Float64
    xmax_active::Float64
    ymin_active::Float64
    ymax_active::Float64
end

function PlotView(plot_dict::PlotDictType)
    params = plot_dict["general_parameters"]
    default_tuple = (0.0, 0.0, 0.0, 0.0)
    return PlotView(
        use_close_up=get(params, "use_close_up", 0),
        xmin_normal=get(params, "dimensions", default_tuple)[1],
        xmax_normal=get(params, "dimensions", default_tuple)[2],
        ymin_normal=get(params, "dimensions", default_tuple)[3],
        ymax_normal=get(params, "dimensions", default_tuple)[4],
        xmin_close_up=get(params, "dim_close_up", default_tuple)[1],
        xmax_close_up=get(params, "dim_close_up", default_tuple)[2],
        ymin_close_up=get(params, "dim_close_up", default_tuple)[3],
        ymax_close_up=get(params, "dim_close_up", default_tuple)[4],
        xmin_active=get(params, "dimensions", default_tuple)[1],
        xmax_active=get(params, "dimensions", default_tuple)[2],
        ymin_active=get(params, "dimensions", default_tuple)[3],
        ymax_active=get(params, "dimensions", default_tuple)[4],
    )
end

function update_active_zoom!(view::PlotView)::Nothing
    if view.use_close_up == 0
        view.xmin_active = view.xmin_normal
        view.xmax_active = view.xmax_normal
        view.ymin_active = view.ymin_normal
        view.ymax_active = view.ymax_normal
    else
        view.xmin_active = view.xmin_close_up
        view.xmax_active = view.xmax_close_up
        view.ymin_active = view.ymin_close_up
        view.ymax_active = view.ymax_close_up
    end
    return nothing
end

function get_active_dimensions(view::PlotView)::Tuple{Float64, Float64, Float64, Float64}
    return (
        view.xmin_active, view.xmax_active,
        view.ymin_active, view.ymax_active
    )
end

end # module