module ModelStructureManager

include("utils/SmoothSurface.jl")
include("utils/TopAndBottom.jl")
include("utils/TopoFromMarkers.jl")
include("utils/LayerAverage.jl")
include("utils/AverageDensity.jl")
include("utils/LocalIsostasy.jl")

include("containers/InterfaceDepthsContainer.jl")
include("containers/LayerThicknessContainer.jl")
include("containers/LayerMaterialIDsContainer.jl")
include("containers/MarkerDataContainer.jl")
include("containers/GridDataContainer.jl")

using .InterfaceDepthsContainer
using .LayerThicknessContainer
using .LayerMaterialIDsContainer
using .MarkerDataContainer
using .GridDataContainer


mutable struct ModelStructure
    marker_data::MarkerDataContainer.MarkerData
    grid_data::GridDataContainer.GridData
    layer_matids::LayerMaterialIDsContainer.LayerMaterialIDs
    interface_depths::InterfaceDepthsContainer.InterfaceDepths
    layer_thickness::LayerThicknessContainer.LayerThickness
    search_radius::Float64
end

function ModelStructure(
    marknum::Int,
    marker_matid::Vector{Int16},
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    xnum::Int,
    gridx_b::Vector{Float64},
    xstp_avg::Float64,
    ysize::Float64,
    matids_sticky::Vector{Int16},
    matids_sticky_water::Vector{Int16},
    matids_rock::Vector{Int16},
    matids_crust::Vector{Int16},
    matids_lithosphere::Vector{Int16},
    matids_asthenosphere::Vector{Int16},
    matids_sediments::Union{Vector{Int16}, Nothing}=nothing
)
    marker_data = MarkerDataContainer.MarkerData(
        marknum=marknum,
        marker_matid=marker_matid,
        marker_x=marker_x,
        marker_y=marker_y
    )
    
    grid_data = GridDataContainer.GridData(
        xnum=xnum,
        gridx_b=gridx_b,
        xstp_avg=xstp_avg,
        ysize=ysize
    )

    search_radius = grid_data.xstp_avg / 2.0

    layer_matids = LayerMaterialIDsContainer.LayerMaterialIDs()
    LayerMaterialIDsContainer.set_material_ids!(
        layer_matids,
        matids_sticky,
        matids_sticky_water,
        matids_rock,
        matids_crust,
        matids_lithosphere,
        matids_asthenosphere,
        matids_sediments
    )

    interface_depths = InterfaceDepthsContainer.InterfaceDepths()
    layer_thickness = LayerThicknessContainer.LayerThickness()

    return ModelStructure(
        marker_data,
        grid_data,
        layer_matids,
        interface_depths,
        layer_thickness,
        search_radius
    )
end

function calculate_structure!(data::ModelStructure)
    calculate_interfaces!(data)
    calculate_layer_thicknesses!(data)
end

function calculate_interfaces!(data::ModelStructure)
    calc_top_and_base_of_sticky_layer(data)
    calc_top_and_base_of_rock_layer(data)
    InterfaceDepthsContainer.calc_ytopo!(data.interface_depths, data.grid_data.xnum)
    InterfaceDepthsContainer.calc_ytopo_smooth!(data.interface_depths, 10, data.grid_data.xnum)
    calc_top_and_base_of_crust_layer(data)
    calc_top_and_base_of_mantle_lithosphere_layer(data)
    InterfaceDepthsContainer.calc_ymoho!(data.interface_depths, data.grid_data.xnum)
    calc_top_and_base_of_asthenosphere_layer(data)
    InterfaceDepthsContainer.calc_ylith_base!(data.interface_depths, data.grid_data.xnum)
    calc_top_and_base_of_sediment_layer(data)
end

function calculate_layer_thicknesses!(data::ModelStructure)
    data.layer_thickness.thick_sticky = calc_sticky_thickness(
        data.grid_data.xnum, data.interface_depths.ytopo)
    data.layer_thickness.thick_crust = calc_crust_thickness(
        data.grid_data.xnum, data.interface_depths.ytopo,
        data.interface_depths.ymoho)
    data.layer_thickness.thick_mantle_lith = calc_mantle_lith_thickness(
        data.grid_data.xnum, data.interface_depths.ymoho,
        data.interface_depths.ylith_base)
    data.layer_thickness.thick_asthenosphere = calc_asthenosphere_thickness(
        data.grid_data.xnum, data.interface_depths.ylith_base,
        data.grid_data.ysize)
    data.layer_thickness.thick_sediments = (
        data.interface_depths.ybottoms_sediments -
        data.interface_depths.ytops_sediments)
end

function get_basic_structure(
    data::ModelStructure
)::Tuple{
    Vector{Float64}, # ytopo
    Vector{Float64}, # ymoho
    Vector{Float64}, # ytopo_smooth
    Vector{Float64}, # ylith_base
    Vector{Float64}, # thick_sticky
    Vector{Float64}, # thick_crust
    Vector{Float64}, # thick_mantle_lith
    Vector{Float64}, # thick_asthenosphere
    Vector{Float64}, # thick_sediments
}
    return (
        data.interface_depths.ytopo,
        data.interface_depths.ymoho,
        data.interface_depths.ytopo_smooth,
        data.interface_depths.ylith_base,
        data.layer_thickness.thick_sticky,
        data.layer_thickness.thick_crust,
        data.layer_thickness.thick_mantle_lith,
        data.layer_thickness.thick_asthenosphere,
        data.layer_thickness.thick_sediments
    )
end

function calculate_water_depth(data::ModelStructure)
    calc_top_and_base_of_sticky_air_layer(data)
    calc_top_and_base_of_sticky_water_layer(data)
    calc_top_and_base_of_rock_layer(data)
    water_depth = InterfaceDepthsContainer.calculate_water_depth(data.interface_depths)
    return water_depth
end

function calc_top_and_base_of_sediment_layer(data::ModelStructure)
    (data.interface_depths.ytops_sediments,
     data.interface_depths.ybottoms_sediments) = 
        get_top_and_bottom(data, data.layer_matids.matids_sediments)
end

function calc_top_and_base_of_sticky_air_layer(data::ModelStructure)
    (data.interface_depths.ytops_sticky_air,
     data.interface_depths.ybottoms_sticky_air) = 
        get_top_and_bottom(data, data.layer_matids.matids_sticky)
end

function calc_top_and_base_of_sticky_water_layer(data::ModelStructure)
    (data.interface_depths.ytops_sticky_water,
     data.interface_depths.ybottoms_sticky_water) = 
        get_top_and_bottom(data, data.layer_matids.matids_sticky_water)
end

function calc_top_and_base_of_sticky_layer(data::ModelStructure)
    (data.interface_depths.ytops_sticky,
     data.interface_depths.ybottoms_sticky) = 
        get_top_and_bottom(data, data.layer_matids.matids_sticky)
end

function calc_top_and_base_of_rock_layer(data::ModelStructure)
    (data.interface_depths.ytops_rock,
     data.interface_depths.ybottoms_rock) = 
        get_top_and_bottom(data, data.layer_matids.matids_rock)
end

function calc_top_and_base_of_crust_layer(data::ModelStructure)
    (data.interface_depths.ytops_crust,
     data.interface_depths.ybottoms_crust) = 
        get_top_and_bottom(data, data.layer_matids.matids_crust)
    InterfaceDepthsContainer.clean_crustal_depths!(
        data.interface_depths, data.grid_data.xnum)
end

function calc_top_and_base_of_mantle_lithosphere_layer(data::ModelStructure)
    (data.interface_depths.ytops_mantle_lith,
     data.interface_depths.ybottoms_mantle_lith) = 
        get_top_and_bottom(data, data.layer_matids.matids_lithosphere)
    InterfaceDepthsContainer.clean_lithosphere_depths!(
        data.interface_depths, data.grid_data.xnum)
end

function calc_top_and_base_of_asthenosphere_layer(data::ModelStructure)
    (data.interface_depths.ytops_asthenosphere,
     data.interface_depths.ybottoms_asthenosphere) = 
        get_top_and_bottom(data, data.layer_matids.matids_asthenosphere)
    InterfaceDepthsContainer.clean_asthenosphere_depths!(
        data.interface_depths, data.grid_data.xnum)
end

function get_top_and_bottom(
    data::ModelStructure,
    matids_layer::Vector{Int16}
)
    (ytops, ybottoms) = TopAndBottom.calculate_top_and_bottom_of_layer_opt(
        matids_layer,
        data.marker_data.marker_matid,
        data.marker_data.marker_x,
        data.marker_data.marker_y,
        data.grid_data.gridx_b,
        data.search_radius
    )
    return ytops, ybottoms
end

function calc_sticky_thickness(
    xnum::Int,
    ytopo::Vector{Float64}
)::Vector{Float64}
    thick_sticky = zeros(Float64, xnum)
    for i in 1:xnum
        thick_sticky[i] = ytopo[i]
    end
    return thick_sticky
end

function calc_crust_thickness(
    xnum::Int,
    ytopo::Vector{Float64},
    ymoho::Vector{Float64}
)::Vector{Float64}
    thick_crust = zeros(Float64, xnum)
    for i in 1:xnum
        thick_crust[i] = ymoho[i] - ytopo[i]
    end
    return thick_crust
end

function calc_mantle_lith_thickness(
    xnum::Int,
    ymoho::Vector{Float64},
    ylith_base::Vector{Float64}
)::Vector{Float64}
    thick_mantle_lith = zeros(Float64, xnum)
    for i in 1:xnum
        thick_mantle_lith[i] = ylith_base[i] - ymoho[i]
    end
    return thick_mantle_lith
end

function calc_asthenosphere_thickness(
    xnum::Int,
    ylith_base::Vector{Float64},
    ysize::Float64
)::Vector{Float64}
    thick_asthenosphere = zeros(Float64, xnum)
    for i in 1:xnum
        thick_asthenosphere[i] = ysize - ylith_base[i]
    end
    return thick_asthenosphere
end

end # module 