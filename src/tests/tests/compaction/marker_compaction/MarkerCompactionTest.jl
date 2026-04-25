module MarkerCompactionTest

import Plots
import EarthBox.Compaction: CompactionCorrection
import EarthBox.Compaction.MarkerCompaction: calculate_marker_swarm_indices
import EarthBox.Compaction.MarkerCompaction: calculate_x_sorted_swarm_indices_from_marker_x
import EarthBox.Compaction.MarkerCompaction: compact_sediment_and_advect_markers
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_layer_opt
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_swarm
import EarthBox.ModelStructureManager.TopAndBottom: calculate_top_and_bottom_of_swarm_opt
import EarthBox.Markers.MarkerCoordinates.InitManager.MarkersRandomized: calc_coordinate_random

function run_test()
    figsize = (10, 5)
    make_plots = true
    ymin_plot = 8000.0
    ymax_plot = 18_000.0
    marker_size = 0.5

    x_domain = 500_000.0
    println("x_domain: ", x_domain)
    y_domain = 18_000.0 # 18_000.0
    println("y_domain: ", y_domain)
    xmin_sed = x_domain * 0.25
    xmax_sed = x_domain * 0.75

    dx_marker = 50.0
    dy_marker = 50.0

    sticky_thickness = 10_000.0
    water_thickness = 1_000.0
    sediment_thickness = 4000.0

    dx_topo = 250.0
    search_radius = dx_topo / 2.0
    nsmooth_top_bottom = 2

    new_sediment_thickness = 200.0
    porosity_initial_transport = 0.5
    decay_depth_transport = 2000.0

    t1 = time()
    (
        markers_x, markers_y, markers_matid,
        markers_porosity_initial, marker_decay_depth, markers_max_burial_depth,
        topo_gridx, topo_gridy, new_sediment_thickness_gridx,
        markers_indices_sedimentary_basin, markers_indices_sticky,
        top_sed, bottom_sed, sediment_thickness_gridx, sticky_thickness_gridx,
        x_sorted_marker_indices_sedimentary_basin
    ) = setup_problem(
        x_domain, y_domain, xmin_sed, xmax_sed, dx_marker, dy_marker,
        sticky_thickness, water_thickness, sediment_thickness, dx_topo,
        new_sediment_thickness, search_radius
    )

    t2 = time()
    println("Time to setup problem: $(t2-t1) seconds")

    if make_plots
        t1 = time()
        plot_markers(
            markers_x, markers_y, markers_matid,
            topo_gridx, topo_gridy, top_sed, bottom_sed,
            name="markers_before_compaction.png", marker_size=marker_size,
            ymin_plot=ymin_plot, ymax_plot=ymax_plot, figsize=figsize
        )
        t2 = time()
        println("Time to plot markers before compaction: $(t2-t1) seconds")
    end

    test_type = "correction" # "compaction" or "correction"
    if test_type == "compaction"
        t1 = time()
        compact_sediment_and_advect_markers(
            topo_gridx,
            topo_gridy,
            sediment_thickness_gridx,
            sticky_thickness_gridx,
            new_sediment_thickness_gridx,
            Float32.(markers_porosity_initial),
            Float32.(marker_decay_depth),
            Float32.(markers_max_burial_depth),
            markers_indices_sedimentary_basin,
            markers_indices_sticky,
            x_sorted_marker_indices_sedimentary_basin,
            markers_x,
            markers_y,
            search_radius,
            nsmooth_top_bottom
        )
        t2 = time()
        println("Time to compact sediment markers: $(t2-t1) seconds")
        t1 = time()
        if make_plots
            plot_markers(
                markers_x, markers_y, markers_matid,
                topo_gridx, topo_gridy, top_sed, bottom_sed,
                name="markers_after_compaction.png", marker_size=marker_size,
                ymin_plot=ymin_plot, ymax_plot=ymax_plot, figsize=figsize
            )
        end
        t2 = time()
        println("Time to plot markers after compaction: $(t2-t1) seconds")

    elseif test_type == "correction"
        topo_gridy_initial = copy(topo_gridy)
        sediment_thickness_initial = copy(sediment_thickness_gridx)
        sticky_thickness_initial = copy(sticky_thickness_gridx)
        topo_gridy_transport = topo_gridy .- new_sediment_thickness_gridx

        nsteps = 3
        for istep in 1:nsteps
            println("Working on correction step: $istep")

            println(">> Number of sticky markers: $(length(markers_indices_sticky))")
            println(">> Number of sediment markers: $(length(markers_indices_sedimentary_basin))")
            println(">> Total number of markers: $(length(markers_x))")
            println(">> Number of topography nodes: $(length(topo_gridx))")

            t1 = time()
            (
                topo_gridy_corrected,
                total_sediment_thickness_corrected,
                new_thickness_decompacted
            ) = CompactionCorrection.decompact_new_sediment_and_compact_markers(
                porosity_initial_transport, decay_depth_transport,
                topo_gridx, topo_gridy_initial, topo_gridy_transport,
                markers_x, markers_y,
                Float32.(markers_porosity_initial), Float32.(marker_decay_depth), Float32.(markers_max_burial_depth), #
                sediment_thickness_initial, sticky_thickness_initial,
                markers_indices_sedimentary_basin, markers_indices_sticky,
                x_sorted_marker_indices_sedimentary_basin,
                search_radius
            )
            println("Min topo_gridy_corrected: $(minimum(topo_gridy_corrected))")
            println("Max topo_gridy_corrected: $(maximum(topo_gridy_corrected))")
            println("Length of topo_gridy_corrected: $(length(topo_gridy_corrected))")

            println("Min thickness decompacted: $(minimum(new_thickness_decompacted))")
            println("Max thickness decompacted: $(maximum(new_thickness_decompacted))")
            t2 = time()
            println("Time to apply compaction correction using markers: $(t2-t1) seconds")

            t1 = time()
            if make_plots
                plot_correction_step(
                    markers_x, markers_y, markers_matid,
                    topo_gridx, topo_gridy_corrected,
                    istep, marker_size=marker_size,
                    ymin_plot=ymin_plot, ymax_plot=ymax_plot, figsize=figsize
                )
            end
            t2 = time()
            println("Time to plot correction step: $(t2-t1) seconds")

            t1 = time()
            topo_gridy_initial = copy(topo_gridy_corrected)
            topo_gridy_transport = topo_gridy_initial .- new_sediment_thickness_gridx
            sediment_thickness_initial = copy(total_sediment_thickness_corrected)
            sticky_thickness_initial = calculate_sticky_thickness(
                markers_matid, markers_x, markers_y, topo_gridx)
            t2 = time()
            println("Time to update variables: $(t2-t1) seconds")
        end
    end
end

""" Setup the problem for compaction of sedimentary basin. """
function setup_problem(
    x_domain::Float64,
    y_domain::Float64,
    xmin_sed::Float64,
    xmax_sed::Float64,
    dx_marker::Float64,
    dy_marker::Float64,
    sticky_thickness::Float64,
    water_thickness::Float64,
    sediment_thickness::Float64,
    dx_topo::Float64,
    new_sediment_thickness::Float64,
    search_radius::Float64
)::Tuple{
    Vector{Float64}, Vector{Float64}, Vector{Int16},
    Vector{Float32}, Vector{Float32}, Vector{Float32},
    Vector{Float64}, Vector{Float64}, Vector{Float64},
    Vector{Int64}, Vector{Int64},
    Vector{Float64}, Vector{Float64},
    Vector{Float64}, Vector{Float64},
    Vector{Int64}
}
    t1 = time()
    (
        markers_x, markers_y, markers_matid,
        markers_porosity_initial, marker_decay_depth, markers_max_burial_depth
    ) = make_markers(
        x_domain, y_domain, xmin_sed, xmax_sed, dx_marker, dy_marker,
        sticky_thickness, water_thickness, sediment_thickness
    )

    topo_gridx = collect(0.0:dx_topo:x_domain-dx_topo)
    println("topo_gridx range: [$(minimum(topo_gridx)), $(maximum(topo_gridx))]")
    topo_gridy = zeros(length(topo_gridx))
    new_sediment_thickness_gridx = define_new_sediment_thickness(
        topo_gridx, xmin_sed, xmax_sed, new_sediment_thickness)

    define_basin_in_topo_grid(
        topo_gridx, topo_gridy, sticky_thickness,
        water_thickness, xmin_sed, xmax_sed
    )

    sedimentary_basin_ids = get_sedimentary_basin_ids()
    sticky_ids = get_sticky_ids()
    t2 = time()
    println("Time to setup markers: $(t2-t1) seconds")

    t1 = time()
    markers_indices_sedimentary_basin = calculate_marker_swarm_indices(
        markers_x, markers_matid, sedimentary_basin_ids
    )

    markers_indices_sticky = calculate_marker_swarm_indices(
        markers_x, markers_matid, sticky_ids
    )
    t2 = time()
    println("Time to calculate marker swarm indices: $(t2-t1) seconds")

    t1 = time()
    x_sorted_marker_indices_sticky = calculate_x_sorted_swarm_indices_from_marker_x(
        markers_x, markers_indices_sticky)
    t2 = time()
    println(">> Time to x-sort sticky swarm indices (s): $(t2-t1)")

    t1 = time()
    x_sorted_marker_indices_sedimentary_basin = calculate_x_sorted_swarm_indices_from_marker_x(
        markers_x, markers_indices_sedimentary_basin)
    t2 = time()
    println(">> Time to x-sort sedimentary basin swarm indices (s): $(t2-t1)")

    t1 = time()
    top_sed, bottom_sed = calculate_top_and_bottom_of_layer_opt(
        sedimentary_basin_ids, markers_matid,
        markers_x, markers_y, topo_gridx, search_radius
    )
    t2 = time()
    println("Time to calculate top and bottom of sedimentary basin: $(t2-t1) seconds")

    sed_thickness_gridx = bottom_sed .- top_sed

    t1 = time()
    sticky_thickness_gridx = calculate_sticky_thickness(
        markers_matid, markers_x, markers_y, topo_gridx
    )
    t2 = time()
    println("Time to calculate sticky thickness: $(t2-t1) seconds")

    sticky_thickness_gridx = zeros(length(topo_gridx))
    use_swarm_approach = true
    if use_swarm_approach
        t1 = time()
        sticky_thickness_gridx = calculate_sticky_thickness_using_swarm(
            markers_x, markers_y, topo_gridx, markers_indices_sticky,
            search_radius
        )
        t2 = time()
        println("Time to calculate sticky thickness (swarm): $(t2-t1) seconds")

        t1 = time()
        sticky_thickness_gridx = calculate_sticky_thickness_using_swarm_opt(
            markers_x, markers_y, topo_gridx,
            x_sorted_marker_indices_sticky, search_radius
        )
        t2 = time()
        println("Time to calculate sticky thickness (swarm opt): $(t2-t1) seconds")
    end

    return (
        markers_x, markers_y, markers_matid,
        markers_porosity_initial, marker_decay_depth, markers_max_burial_depth,
        topo_gridx, topo_gridy, new_sediment_thickness_gridx,
        markers_indices_sedimentary_basin, markers_indices_sticky,
        top_sed, bottom_sed,
        sed_thickness_gridx, sticky_thickness_gridx,
        x_sorted_marker_indices_sedimentary_basin
    )
end

function calculate_sticky_thickness(
    markers_matid::Vector{Int16},
    markers_x::Vector{Float64},
    markers_y::Vector{Float64},
    topo_gridx::Vector{Float64}
)::Vector{Float64}
    sticky_ids = get_sticky_ids()
    search_radius = topo_gridx[2] - topo_gridx[1]
    top_sticky, bottom_sticky = calculate_top_and_bottom_of_layer_opt(
        sticky_ids, markers_matid,
        markers_x, markers_y, topo_gridx, search_radius
    )
    sticky_thickness = bottom_sticky - top_sticky
    return sticky_thickness
end

function calculate_sticky_thickness_using_swarm(
    markers_x::Vector{Float64},
    markers_y::Vector{Float64},
    topo_gridx::Vector{Float64},
    markers_indices_sticky::Vector{Int64},
    search_radius::Float64
)::Vector{Float64}
    top_sticky, bottom_sticky = calculate_top_and_bottom_of_swarm(
        markers_indices_sticky,
        markers_x, markers_y, topo_gridx, search_radius
    )
    sticky_thickness = bottom_sticky - top_sticky
    return sticky_thickness
end

function calculate_sticky_thickness_using_swarm_opt(
    markers_x::Vector{Float64},
    markers_y::Vector{Float64},
    topo_gridx::Vector{Float64},
    x_sorted_marker_indices_sticky::Vector{Int64},
    search_radius::Float64
)::Vector{Float64}
    top_sticky, bottom_sticky = calculate_top_and_bottom_of_swarm_opt(
        x_sorted_marker_indices_sticky,
        markers_x, markers_y, topo_gridx, search_radius
    )
    sticky_thickness = bottom_sticky - top_sticky
    return sticky_thickness
end

function get_matids()::Tuple{Int16, Int16, Int16, Int16, Int16}
    matid_sticky = Int16(1)
    matid_water = Int16(2)
    matid_compacting_sediment = Int16(3)
    matid_middle_sediment_layers = Int16(4)
    matid_crust = Int16(5)
    return (
        matid_sticky, matid_water, matid_compacting_sediment,
        matid_middle_sediment_layers, matid_crust
    )
end

function get_sedimentary_basin_ids()::Vector{Int16}
    (
        _, _, matid_compacting_sediment,
        matid_middle_sediment_layers, _
    ) = get_matids()
    sedimentary_basin_matids = [
        matid_compacting_sediment, matid_middle_sediment_layers
    ]
    return sedimentary_basin_matids
end

function get_sticky_ids()::Vector{Int16}
    (
        matid_sticky, matid_water, _, _, _
    ) = get_matids()
    sticky_ids = [matid_sticky, matid_water]
    return sticky_ids
end

function make_markers(
    x_domain::Float64,
    y_domain::Float64,
    xmin_sed::Float64,
    xmax_sed::Float64,
    dx_marker::Float64,
    dy_marker::Float64,
    sticky_thickness::Float64,
    water_thickness::Float64,
    sediment_thickness::Float64,
    porosity_initial::Float64 = 0.5,
    decay_depth::Float64 = 2000.0
)::Tuple{
    Vector{Float64}, Vector{Float64}, Vector{Int16},
    Vector{Float32}, Vector{Float32}, Vector{Float32}
}
    (
        matid_sticky, matid_water, matid_compacting_sediment,
        matid_middle_sediment_layers, matid_crust
    ) = get_matids()

    nx = floor(Int, x_domain/dx_marker) + 1
    ny = floor(Int, y_domain/dy_marker) + 1
    nmarkers = nx * ny

    ymin_sed = sticky_thickness + water_thickness
    ymax_sed = ymin_sed + sediment_thickness

    ymin_middle_sed_layer = ymin_sed + 0.25 * sediment_thickness
    ymax_middle_sed_layer = ymin_sed + 0.50 * sediment_thickness

    use_random = true
    if !use_random
        markers_x, markers_y = define_marker_coordinates(nx, ny, dx_marker, dy_marker)
    else
        markers_x, markers_y = define_marker_coordinates_random(nx, ny, dx_marker, dy_marker)
    end

    markers_matid = fill(matid_crust, nmarkers)
    markers_porosity_initial = fill(Float32(porosity_initial), nmarkers)
    marker_decay_depth = fill(Float32(decay_depth), nmarkers)

    define_background_crustal_material_ids(
        markers_matid, markers_y, sticky_thickness, matid_sticky, matid_crust
    )

    define_material_ids_of_water_layer(
        markers_matid, markers_x, markers_y,
        matid_water, sticky_thickness,
        xmin_sed, xmax_sed, ymin_sed
    )

    define_material_ids_of_sedimentary_basin(
        markers_matid, markers_x, markers_y,
        matid_compacting_sediment,
        xmin_sed, xmax_sed, ymin_sed, ymax_sed
    )

    define_material_ids_of_middle_layer(
        markers_matid,
        markers_porosity_initial,
        marker_decay_depth,
        markers_y,
        matid_compacting_sediment,
        matid_middle_sediment_layers,
        ymin_middle_sed_layer,
        ymax_middle_sed_layer,
        0.0, # Set to non-compacting material
        1.0
    )

    markers_max_burial_depth = zeros(Float32, nmarkers)
    return (
        markers_x, markers_y, markers_matid, markers_porosity_initial,
        marker_decay_depth, markers_max_burial_depth
    )
end

""" Define the sedimentary basin in the topography grid. """
function define_basin_in_topo_grid(
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    sticky_thickness::Float64,
    water_thickness::Float64,
    xmin_sed::Float64,
    xmax_sed::Float64
)::Nothing
    ntopo = length(topo_gridx)
    for i in 1:ntopo
        topo_gridy[i] = sticky_thickness
        if xmin_sed <= topo_gridx[i] <= xmax_sed
            topo_gridy[i] = sticky_thickness + water_thickness
        end
    end
    return nothing
end

function define_new_sediment_thickness(
    topo_gridx::Vector{Float64},
    xmin_sed::Float64,
    xmax_sed::Float64,
    new_sediment_thickness::Float64
)::Vector{Float64}
    new_sediment_thickness_gridx = fill(new_sediment_thickness, length(topo_gridx))
    ntopo = length(topo_gridx)
    for i in 1:ntopo
        xgrid = topo_gridx[i]
        if xgrid < xmin_sed || xgrid >= xmax_sed
            new_sediment_thickness_gridx[i] = 0.0
        end
    end
    return new_sediment_thickness_gridx
end

function define_marker_coordinates(
    nx::Int64,
    ny::Int64,
    dx::Float64,
    dy::Float64
)::Tuple{Vector{Float64}, Vector{Float64}}
    nmarkers = nx * ny
    marker_x = zeros(nmarkers)
    marker_y = zeros(nmarkers)
    # Create marker coordinates using dx and dy
    icount = 1
    for j in 1:ny
        for i in 1:nx
            marker_x[icount] = Float64(i-1) * dx
            marker_y[icount] = Float64(j-1) * dy
            icount += 1
        end
    end
    return marker_x, marker_y
end

function define_marker_coordinates_random(
    nx::Int64,
    ny::Int64,
    dx::Float64,
    dy::Float64
)::Tuple{Vector{Float64}, Vector{Float64}}
    nmarkers = nx * ny
    marker_x = zeros(nmarkers)
    marker_y = zeros(nmarkers)
    random_numbers = rand(nmarkers)
    # Create marker coordinates using dx and dy
    icount = 1
    for j in 1:ny
        for i in 1:nx
            marker_x[icount] = calc_coordinate_random(i-1, dx, random_numbers[icount])
            marker_y[icount] = calc_coordinate_random(j-1, dy, random_numbers[icount])
            icount += 1
        end
    end
    return marker_x, marker_y
end

function define_background_crustal_material_ids(
    markers_matid::Vector{Int16},
    marker_y::Vector{Float64},
    sticky_thickness::Float64,
    matid_sticky::Int16,
    matid_crust::Int16
)::Nothing
    nmarkers = length(markers_matid)
    for i in 1:nmarkers
        if marker_y[i] <= sticky_thickness
            markers_matid[i] = matid_sticky
        else
            markers_matid[i] = matid_crust
        end
    end
    return nothing
end

function define_material_ids_of_sedimentary_basin(
    markers_matid::Vector{Int16},
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    matid_sedimentary_basin::Int16,
    xmin_sed::Float64,
    xmax_sed::Float64,
    ymin_sed::Float64,
    ymax_sed::Float64
)::Nothing
    nmarkers = length(markers_matid)
    for i in 1:nmarkers
        if xmin_sed <= marker_x[i] <= xmax_sed && ymin_sed <= marker_y[i] < ymax_sed
            markers_matid[i] = matid_sedimentary_basin
        end
    end
    return nothing
end

function define_material_ids_of_water_layer(
    markers_matid::Vector{Int16},
    marker_x::Vector{Float64},
    marker_y::Vector{Float64},
    matid_water::Int16,
    sticky_thickness::Float64,
    xmin_sediment::Float64,
    xmax_sediment::Float64,
    ymin_sediment::Float64
)::Nothing
    nmarkers = length(markers_matid)
    for i in 1:nmarkers
        if xmin_sediment <= marker_x[i] <= xmax_sediment
            if sticky_thickness < marker_y[i] < ymin_sediment
                markers_matid[i] = matid_water
            end
        end
    end
    return nothing
end

function define_material_ids_of_middle_layer(
    markers_matid::Vector{Int16},
    marker_porosity_initial::Vector{Float32},
    marker_decay_depth::Vector{Float32},
    marker_y::Vector{Float64},
    matid_sedimentary_basin::Int16,
    matid_middle_layer::Int16,
    y_min::Float64,
    y_max::Float64,
    porosity_initial::Float64,
    decay_depth::Float64
)::Nothing
    nmarkers = length(markers_matid)
    for i in 1:nmarkers
        if y_min < marker_y[i] < y_max
            if markers_matid[i] == matid_sedimentary_basin
                markers_matid[i] = matid_middle_layer
                # Set marker to non-compacting material
                marker_porosity_initial[i] = porosity_initial
                marker_decay_depth[i] = decay_depth
            end
        end
    end
    return nothing
end

function plot_markers(
    markers_x::Vector{Float64},
    markers_y::Vector{Float64},
    markers_matid::Vector{Int16},
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    top_sed::Vector{Float64},
    bottom_sed::Vector{Float64};
    name::String = "markers.png",
    marker_size::Float64 = 10.0,
    plot_top_bottom_sed::Bool = true,
    ymin_plot::Float64 = 0.0,
    ymax_plot::Float64 = 18_000.0,
    figsize::Tuple{Int, Int} = (10, 5)
)::Nothing
    # Plot the markers with different colors for different material ids
    dpi = 150
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)

    # Print min/max values of marker arrays
    println("markers_x range: [$(minimum(markers_x)), $(maximum(markers_x))]")
    println("markers_y range: [$(minimum(markers_y)), $(maximum(markers_y))]") 
    println("markers_matid range: [$(minimum(markers_matid)), $(maximum(markers_matid))]")

    # Make matid colored red
    markercolor = [matid == 3 ? :red : :blue for matid in markers_matid]
    p = Plots.scatter(
        markers_x, markers_y, markercolor=markercolor, markersize=marker_size,
        markerstrokewidth=0.0, markerstrokecolor=:transparent,
        markershape=:rect, size=figsize_pixels, dpi=dpi, xlabel="x", ylabel="y",
        ylims=(ymin_plot, ymax_plot), aspect_ratio=:auto, yflip=true
        )
    Plots.plot!(p, topo_gridx, topo_gridy, color=:red, label="Topo", linewidth=2)
    if plot_top_bottom_sed
        Plots.plot!(p, topo_gridx, top_sed, color=:green, linewidth=2, label="Top Sed")
        Plots.plot!(
            p, topo_gridx, bottom_sed, color=:black, linestyle=:dash, 
            linewidth=2, label="Bottom Sed"
            )
    end
    Plots.savefig(p, name)
    return nothing
end

function plot_correction_step(
    markers_x::Vector{Float64},
    markers_y::Vector{Float64},
    markers_matid::Vector{Int16},
    topo_gridx::Vector{Float64},
    topo_gridy::Vector{Float64},
    istep::Int64;
    marker_size::Float64 = 10.0,
    ymin_plot::Float64 = 0.0,
    ymax_plot::Float64 = 18_000.0,
    figsize::Tuple{Int, Int} = (10, 5)
)::Nothing
    # Plot the markers with different colors for different material ids
    dpi = 150
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)
    # Make matid colored red
    markercolor = [matid == 3 ? :red : :blue for matid in markers_matid]
    p = Plots.scatter(
        markers_x, markers_y, markercolor=markercolor, markersize=marker_size,
        markerstrokewidth=0.0, markerstrokecolor=:transparent,
        markershape=:rect, size=figsize_pixels, dpi=dpi, xlabel="x", ylabel="y",
        ylims=(ymin_plot, ymax_plot), aspect_ratio=:auto, yflip=true
        )
    Plots.plot!(p, topo_gridx, topo_gridy, color=:red, linewidth=2)
    Plots.savefig(p, "compaction_correction_step_$(istep).png")
    return nothing
end

function plot_sediment_thickness_corrected(
    topo_gridx::Vector{Float64},
    sediment_thickness_corrected::Vector{Float64},
    istep::Int64
)::Nothing
    dpi = 150
    figsize = (5, 5)
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)
    p = Plots.plot(
        xlabel="x",
        ylabel="Sediment Thickness",
        ylims=(0, 5000.0),
        size=figsize_pixels
    )
    Plots.plot!(p, topo_gridx, sediment_thickness_corrected)
    Plots.savefig(p, "sediment_thickness_corrected_$(istep).png")
    return nothing
end

function plot_sticky_thickness(
    topo_gridx::Vector{Float64},
    sticky_thickness::Vector{Float64},
    istep::Int64
)::Nothing
    dpi = 150
    figsize = (5, 5)
    figsize_pixels = (figsize[1] * dpi, figsize[2] * dpi)
    p = Plots.plot(
        xlabel="x",
        ylabel="Sticky Thickness",
        ylims=(0, 4000.0),
        size=figsize_pixels
    )
    Plots.plot!(p, topo_gridx, sticky_thickness)
    Plots.savefig(p, "sticky_thickness$(istep).png")
    return nothing
end

end # module 