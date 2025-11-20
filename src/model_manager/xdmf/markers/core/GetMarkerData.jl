module GetMarkerData

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer.OutputStandard: OutputLists
import EarthBox.Arrays.ArrayTypes.MarkerArrayFloat1D: MarkerArrayFloat1DState
import EarthBox.ConfigurationManager.OutputConfig: OutputConfigState
import EarthBox.SurfaceProcesses.Sealevel.UpdateSealevel: get_time_dependent_base_level_shift
import ...OutputDTypes: ScalarField
import ...XdmfUtils: get_xdmf_number_type_for_array
import ...XdmfUtils: getoutform

function get_marker_data(
    model::ModelData,
    output_lists::OutputLists,
    output_config::OutputConfigState
)::Dict{String, Any}
    noutput = model.timestep.parameters.output_steps.noutput.value
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    sec_per_Myr = model.conversion.parameters.sec_per_Myr.value

    model_time = timesum/sec_per_Myr
    time_units = "Myr"

    y_sealevel = model.topography.parameters.sealevel.y_sealevel.value
    base_level_shift = get_time_dependent_base_level_shift(model)
    array_object_list = output_lists.marker_array_obj_list

    scalars_on_markers = []

    for obj in array_object_list
        create_output_for_object = output_config.marker_output[obj.name]
        if create_output_for_object
            scalar_field = ScalarField(
                obj.outform.header,
                obj.name,
                get_xdmf_number_type_for_array(obj.array),
                obj.outform.units,
                obj.outform.header,
                "markers",
                getoutform(obj)
            )
            push!(scalars_on_markers, scalar_field)
        end
    end

    jld_markerfile = "markers_$(lpad(noutput, 5, "0")).jld"

    nmarkers, marker_xy_km, marker_ids = make_marker_arrays_for_xdmf_format(
        array_object_list[1],
        array_object_list[2]
    )

    marker_data = Dict{String, Any}(
        "nmarkers" => nmarkers,
        "jld_markerfile" => jld_markerfile,
        "jld_dataname_xy" => "marker_xy_km",
        "jld_dataname_ids" => "marker_ids",
        "marker_id_array" => marker_ids,
        "marker_xy_km_array" => marker_xy_km,
        "scalars_on_markers" => scalars_on_markers,
        "time" => model_time,
        "time_units" => time_units,
        "noutput" => noutput,
        "y_sealevel" => y_sealevel,
        "base_level_shift" => base_level_shift
    )

    return marker_data
end

function make_marker_arrays_for_xdmf_format(
    marker_x::MarkerArrayFloat1DState{Float64},
    marker_y::MarkerArrayFloat1DState{Float64}
)::Tuple{Int, Array{Float64, 2}, Array{Float64, 1}}
    marker_x_m = getoutform(marker_x)
    marker_y_m = -getoutform(marker_y)
    nmarkers = size(marker_x_m, 1)
    marker_xy_km = zeros(Float64, nmarkers, 2)
    marker_ids = zeros(Float64, nmarkers)
    update_marker_xy_km_array(nmarkers, marker_x_m, marker_y_m, marker_xy_km)
    marker_ids = collect(0:nmarkers-1)
    return nmarkers, marker_xy_km, marker_ids
end

function update_marker_xy_km_array(
    nmarkers::Int,
    marker_x_m::Array{Float64, 1},
    marker_y_m::Array{Float64, 1},
    marker_xy_km::Array{Float64, 2}
)::Nothing
    for i in 1:nmarkers
        marker_xy_km[i, 1] = marker_x_m[i]/1000.0
        marker_xy_km[i, 2] = marker_y_m[i]/1000.0
    end
    return nothing
end

end # module 