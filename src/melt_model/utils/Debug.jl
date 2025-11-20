module Debug

import EarthBox.ModelDataContainer: ModelData
import EarthBox.ConversionFuncs: kelvin_to_celsius, celsius_to_kelvin, seconds_to_years
import ..MeltFraction: get_melting_model_parameters

function print_melting_info(model::ModelData)::Nothing
    box_xmid = (243.74 + 256.25)/2.0*1000.0 # meters
    box_ymid = (16.25 + 17.5)/2.0*1000.0 # meters
    box_width = 2500.0 # meters
    box_height = 2500.0 # meters
    imarker_target = -1
    time_to_print_myr = 0.8
    print_max = 10
    temperature_min = celsius_to_kelvin(1200.0)
    timesum = model.timestep.parameters.main_time_loop.timesum.value
    timesum_yr = seconds_to_years(timesum)
    timesum_myr = timesum_yr*1.0e-6

    if timesum_myr >= time_to_print_myr
        icount = 0
        marknum = model.markers.parameters.distribution.marknum.value
        for imarker in 1:marknum
            if icount < print_max
                icount = print_marker_melting_info(
                    model, imarker, icount, timesum_myr,
                    box_xmid, box_ymid, box_width, box_height,
                    temperature_min, imarker_target
                )
            end
        end
    end
    return nothing
end

function print_marker_melting_info(
    model::ModelData,
    imarker::Int64,
    icount::Int64,
    timesum_myr::Float64,
    box_xmid::Float64,
    box_ymid::Float64,
    box_width::Float64,
    box_height::Float64,
    temperature_min::Float64,
    imarker_target::Int64 = -1
)::Int64
    xmin = box_xmid - 0.5*box_width
    xmax = box_xmid + 0.5*box_width
    ymin = box_ymid - 0.5*box_height
    ymax = box_ymid + 0.5*box_height

    # Unpack coordinate arrays
    marker_x = model.markers.arrays.location.marker_x.array
    marker_y = model.markers.arrays.location.marker_y.array
    # Get coordinates
    x_marker = marker_x[imarker]
    y_marker = marker_y[imarker]
    marker_tk = model.markers.arrays.thermal.marker_TK.array
    marker_temp = marker_tk[imarker]

    printit = false
    if imarker_target < 0
        # Use box if imarker_target is negative
        if xmin < x_marker < xmax && ymin < y_marker < ymax
            printit = true
        end
    else
        if imarker == imarker_target
            printit = true
        end
    end

    if marker_temp < temperature_min
        printit = false
    end

    if printit
        # Unpack arrays
        marker_matid = model.markers.arrays.material.marker_matid.array
        marker_meltfrac = model.markers.arrays.melt.marker_meltfrac.array
        marker_extracted_meltfrac = model.markers.arrays.melt.marker_extracted_meltfrac.array
        marker_extractable_meltfrac = model.markers.arrays.melt.marker_extractable_meltfrac.array
        marker_pr = model.markers.arrays.pressure.marker_pr.array
        mat_melting_itypes = model.materials.arrays.mat_melting_itypes.array
        mat_melting = model.materials.arrays.mat_melting.array

        # Get array values
        matid = marker_matid[imarker]
        meltfrac = marker_meltfrac[imarker]
        extracted_meltfrac = marker_extracted_meltfrac[imarker]
        extractable_meltfrac = marker_extractable_meltfrac[imarker]
        marker_pr_gpa = marker_pr[imarker]*1.0e-9

        itype_solidus = mat_melting_itypes[matid, 1]
        itype_liquidus = mat_melting_itypes[matid, 2]
        latent_heat = mat_melting[matid, 1]

        temperature_liquidus, temperature_solidus = get_melting_model_parameters(
            marker_pr_gpa*1e9,
            itype_solidus,
            itype_liquidus
        )

        println(">> Melting info: imarker: ", imarker, " timesum(Myr): ", 
                round(timesum_myr, digits=3))
        println("matid: ", matid, " tempC: ", round(kelvin_to_celsius(marker_temp), digits=2),
                " prGPa: ", round(marker_pr_gpa, digits=2), " x: ", 
                round(x_marker, digits=2), " y: ", round(y_marker, digits=2))
        println("meltfrac: ", round(meltfrac, digits=3), " extracted_meltfrac: ", 
                round(extracted_meltfrac, digits=3), " extractable_meltfrac: ", 
                round(extractable_meltfrac, digits=3))
        println("itype_solidus: ", itype_solidus, " itype_liquidus: ", itype_liquidus,
                " latent_heat: ", round(latent_heat, digits=2), " temp_liquidus(C): ", 
                round(kelvin_to_celsius(temperature_liquidus), digits=2),
                " temp_solidus(C): ", round(kelvin_to_celsius(temperature_solidus), digits=2))
        icount += 1
    end

    return icount
end

end # module 