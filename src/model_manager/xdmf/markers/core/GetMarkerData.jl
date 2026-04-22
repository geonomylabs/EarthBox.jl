module GetMarkerData

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ModelDataContainer.OutputStandard: OutputLists
import EarthBox.ConfigurationManager.OutputConfig: OutputConfigState
import EarthBox.SurfaceProcesses.Sealevel.UpdateSealevel: get_time_dependent_base_level_shift
import ...OutputDTypes: ScalarFieldMeta
import ...XdmfUtils: get_xdmf_number_type_for_array

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

    scalar_metas = ScalarFieldMeta[]
    enabled_marker_objs = []

    for obj in array_object_list
        create_output_for_object = output_config.marker_output[obj.name]
        if create_output_for_object
            push!(scalar_metas, ScalarFieldMeta(
                obj.outform.header,
                obj.name,
                get_xdmf_number_type_for_array(obj.array),
                obj.outform.units,
                obj.outform.header,
                "markers"
            ))
            push!(enabled_marker_objs, obj)
        end
    end

    jld_markerfile = "markers_$(lpad(noutput, 5, "0")).jld"

    nmarkers = length(array_object_list[1].array)

    marker_data = Dict{String, Any}(
        "nmarkers" => nmarkers,
        "marker_x_obj" => array_object_list[1],
        "marker_y_obj" => array_object_list[2],
        "jld_markerfile" => jld_markerfile,
        "jld_dataname_x" => "marker_x_km",
        "jld_dataname_y" => "marker_y_km",
        "jld_dataname_ids" => "marker_ids",
        "scalar_metas" => scalar_metas,
        "enabled_marker_objs" => enabled_marker_objs,
        "time" => model_time,
        "time_units" => time_units,
        "noutput" => noutput,
        "y_sealevel" => y_sealevel,
        "base_level_shift" => base_level_shift
    )

    return marker_data
end

end # module