module GetData

import EarthBox: MathTools
import ...Reader

function get_inputs(
    itime_step::Int,
    main_paths::Dict{String, String}
)::Dict{String, Any}
    file_base_name = "particles"
    (
        marker_x_m,
        marker_y_m,
        marker_matid,
        tmyr
    ) = Reader.read_marker_material_ids(
        main_paths["post_proc_input_path"], file_base_name, itime_step
    )

    file_base_name = "rho_kg_m3"
    (
        gridx_b_km,
        gridy_b_km,
        rho_kg_m3,
        tmyr
    ) = Reader.read_scalar_grid(
        main_paths["post_proc_input_path"], file_base_name, itime_step
    )

    file_base_name = "topo"
    (
        xtopo_marker_chain_full,
        ytopo_marker_chain_full
    ) = Reader.read_topography(
        main_paths["post_proc_input_path"], file_base_name, itime_step
    )

    file_base_name = "topo"
    (
        xtopo_marker_chain_ini_full,
        ytopo_marker_chain_ini_full
    ) = Reader.read_topography(
        main_paths["post_proc_input_path"], file_base_name, 1
    )

    gridx_b = copy(gridx_b_km) .* 1000.0
    gridy_b = copy(gridy_b_km) .* 1000.0
    xnum = size(gridx_b_km, 1)
    ynum = size(gridy_b_km, 1)
    xsize = gridx_b[end]
    ysize = gridy_b[end]
    marknum = length(marker_x_m)
    xstp_avg = xsize / (xnum - 1)

    ytopo_marker_chain = zeros(xnum)
    MathTools.linear_interp_vals!(
        xtopo_marker_chain_full,
        ytopo_marker_chain_full,
        gridx_b,
        ytopo_marker_chain
    )

    ytopo_marker_chain_ini = zeros(xnum)
    MathTools.linear_interp_vals!(
        xtopo_marker_chain_ini_full,
        ytopo_marker_chain_ini_full,
        gridx_b,
        ytopo_marker_chain_ini
    )

    input_dict = Dict{String, Any}(
        "tmyr" => tmyr,
        "gridx_b" => gridx_b,
        "gridy_b" => gridy_b,
        "rho_kg_m3" => rho_kg_m3,
        "xnum" => xnum,
        "ynum" => ynum,
        "xsize" => xsize,
        "ysize" => ysize,
        "xstp_avg" => xstp_avg,
        "marknum" => marknum,
        "marker_x_m" => marker_x_m,
        "marker_y_m" => marker_y_m,
        "marker_matid" => marker_matid,
        "ytopo_marker_chain" => ytopo_marker_chain,
        "ytopo_marker_chain_ini" => ytopo_marker_chain_ini
    )

    return input_dict
end

end # module 