module HotBoxMeltingManager

import EarthBox.EarthBoxDtypes: TestInfoDictType
import EarthBox.EarthBoxDtypes: ModelInputDictType
import EarthBox.EarthBoxDtypes: MaterialsDictType
import ...Reader
import ...TestResults
import ...BenchmarksStruct: Benchmarks

const MATID_SOLIDIFIED_LAYERED_GABBRO = Int16(10)
const MATID_SOLIDIFIED_GABBRO         = Int16(7)
const MATID_SOLIDIFIED_BASALT         = Int16(14)

const TARGET_LAYERED_GABBRO_COUNT = 79721.0
const TARGET_GABBRO_COUNT         = 39224.0
const TARGET_BASALT_COUNT         = 26769.0

const RELATIVE_ERROR_LIMIT_PERCENTAGE = 2.0

function compare_numerical_to_numerical(
    bench::Benchmarks
)::Tuple{Vector{Union{String, Float64}}, String}
    _marker_x_m, _marker_y_m, marker_matid, tmyr = Reader.read_marker_material_ids(
        bench.main_paths["post_proc_input_path"],
        "particles",
        bench.itime_step
    )

    count_layered_gabbro = 0
    count_gabbro = 0
    count_basalt = 0
    for matid in marker_matid
        if matid == MATID_SOLIDIFIED_LAYERED_GABBRO
            count_layered_gabbro += 1
        elseif matid == MATID_SOLIDIFIED_GABBRO
            count_gabbro += 1
        elseif matid == MATID_SOLIDIFIED_BASALT
            count_basalt += 1
        end
    end

    numerical = Float64[
        count_layered_gabbro,
        count_gabbro,
        count_basalt
    ]
    target = Float64[
        TARGET_LAYERED_GABBRO_COUNT,
        TARGET_GABBRO_COUNT,
        TARGET_BASALT_COUNT
    ]

    println("")
    println("Hot Box Melting marker counts at tmyr = ", tmyr)
    println(
        "  Solidified Layered Gabbro (matid 10): ",
        "Numerical = ", count_layered_gabbro,
        " Target = ", Int(TARGET_LAYERED_GABBRO_COUNT)
    )
    println(
        "  Solidified Gabbro          (matid 7):  ",
        "Numerical = ", count_gabbro,
        " Target = ", Int(TARGET_GABBRO_COUNT)
    )
    println(
        "  Solidified Basalt          (matid 14): ",
        "Numerical = ", count_basalt,
        " Target = ", Int(TARGET_BASALT_COUNT)
    )
    println("")

    result, result_msg = TestResults.get_test_results_numerical_vs_analytical(
        numerical,
        target,
        RELATIVE_ERROR_LIMIT_PERCENTAGE / 100.0
    )
    return result, result_msg
end

end # module
