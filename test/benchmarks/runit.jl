import DataStructures: OrderedDict
using EarthBox

function run_benchmarks()
    option_names = BenchmarksManager.option_names
    test_dict = OrderedDict(
        option_names.couette_flow_viscous_heating                          => false,
        option_names.channel_flow_non_newtonian                            => false,
        option_names.channel_flow_variable_conductivity                    => false,
        option_names.channel_flow_non_steady_temperature                   => false,
        option_names.solid_body_rotation                                   => false,
        option_names.rayleigh_taylor_instability                           => false,
        option_names.elastic_slab                                          => false,
        option_names.viscoelastic_stress_buildup                           => false,
        option_names.plasticity_benchmark_kaus10                           => false,
        option_names.box_convection_isoviscous_1a                          => false,
        option_names.viscoelastic_extension                                => false,
        option_names.viscoelastic_extension_asymmetric                     => false,
        option_names.viscoelastic_extension_depth_dependent                => false,
        option_names.viscoelastic_extension_inflow_and_outflow_along_sides => false,
        option_names.viscoelastic_contraction                              => false,
        option_names.viscoelastic_contraction_asymmetric                   => false,
        option_names.simple_sedimentation                                  => false,
        option_names.seafloor_spreading                                    => false,
        option_names.flexure_triangular_hole                               => true,
    )
    mumps_solver_dict = Dict{Symbol, Vector{Any}}(
        option_names.flexure_triangular_hole                               => [true, 8],
        option_names.viscoelastic_extension_inflow_and_outflow_along_sides => [false, 20]
    )
    BenchmarksManager.run_benchmarks(
        test_dict,
        mumps_solver_dict   = mumps_solver_dict,
        run_model           = true,
        run_post_processing = true,
        base_path           = "/mnt/extradrive1",
        old_date_stamp      = nothing,
        make_backup         = false,
        restart_from_backup = false
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_benchmarks()
end 
