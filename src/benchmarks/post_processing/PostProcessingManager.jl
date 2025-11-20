module PostProcessingManager

include("analytical/utils/GetData.jl")
include("analytical/utils/ViscoelasticLithosphericDeformation.jl")
include("analytical/ChannelFlowNonNewtonianManager.jl")
include("analytical/ChannelFlowVariableConductivityManager.jl")
include("analytical/CouetteFlowViscousHeatingManager.jl")
include("analytical/ChannelFlowNonSteadyTemperatureManager.jl")
include("analytical/ViscoElasticStressBuildupManager.jl")
include("analytical/RayleighTaylorInstabilityManager.jl")
include("analytical/SolidBodyRotationManager.jl")
include("analytical/ElasticSlabManager.jl")
include("analytical/ViscoelasticExtension.jl")
include("analytical/ViscoelasticContraction.jl")
include("analytical/flexure_tests/FlexureTriangularHole.jl")
include("numerical/PlasticityBenchmarkK10.jl")
include("analytical/SimpleSedimentationManager.jl")
include("numerical/box_convection/BoxConvectionManager.jl")
include("numerical/seafloor_spreading/SeafloorSpreadingManager.jl")

import .ChannelFlowNonNewtonianManager
import .ChannelFlowVariableConductivityManager
import .CouetteFlowViscousHeatingManager
import .RayleighTaylorInstabilityManager
import .SolidBodyRotationManager
import .ElasticSlabManager
import .SimpleSedimentationManager
import .ViscoelasticExtension
import .ViscoelasticContraction
import ..Options: option_names
import ..BenchmarksStruct: Benchmarks

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.channel_flow_non_newtonian}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ChannelFlowNonNewtonianManager.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.channel_flow_variable_conductivity}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ChannelFlowVariableConductivityManager.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.couette_flow_viscous_heating}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = CouetteFlowViscousHeatingManager.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.channel_flow_non_steady_temperature}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ChannelFlowNonSteadyTemperatureManager.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.viscoelastic_stress_buildup}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ViscoElasticStressBuildupManager.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.rayleigh_taylor_instability}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = RayleighTaylorInstabilityManager.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.solid_body_rotation}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = SolidBodyRotationManager.compare_to_initial_temperature_wave(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.elastic_slab}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ElasticSlabManager.compare_elastic_slab_marker_position_to_initial(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.plasticity_benchmark_kaus10}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = PlasticityBenchmarkK10.compare_numerical_to_numerical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.viscoelastic_extension}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ViscoelasticExtension.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.viscoelastic_extension_asymmetric}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ViscoelasticExtension.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.viscoelastic_extension_depth_dependent}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ViscoelasticExtension.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.viscoelastic_extension_inflow_and_outflow_along_sides}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ViscoelasticExtension.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.viscoelastic_contraction}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ViscoelasticContraction.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.viscoelastic_contraction_asymmetric}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = ViscoelasticContraction.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.flexure_triangular_hole}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = FlexureTriangularHole.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.simple_sedimentation}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = SimpleSedimentationManager.compare_numerical_to_analytical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.box_convection_isoviscous_1a}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = BoxConvectionManager.compare_numerical_to_numerical(bench)
    return result, result_msg
end

function post_proc_func(
    bench::Benchmarks,
    ::Val{option_names.seafloor_spreading}
)::Tuple{Vector{Union{String, Float64}}, String}
    (
        result, result_msg
    ) = SeafloorSpreadingManager.compare_numerical_to_empirical(bench)
    return result, result_msg
end

end # module

