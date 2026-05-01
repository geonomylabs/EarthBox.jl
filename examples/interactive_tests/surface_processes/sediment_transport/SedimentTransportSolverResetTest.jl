module SedimentTransportSolverResetTest

import EarthBox.ConversionFuncs: mm_per_yr_to_meters_per_seconds as mm_yr_to_m_s
import EarthBox.ConversionFuncs: years_to_seconds
import EarthBox.ConversionFuncs: meters_per_year_to_meters_per_seconds as m_yr_to_m_s
import EarthBox.ConversionFuncs: meters_squared_per_year_to_meters_squared_per_second as m2_yr_to_m2_s
import EarthBox.DataStructures: SedimentTransportParameters
import EarthBox.SurfaceProcesses.SedimentTransport: SedimentTransportSolverManager

const _SOLVER_KWARGS = (
    use_collections = true,
    use_print_debug = false,
    use_constant_diffusivity = false,
    use_compaction_correction = true,
    compaction_correction_type = "constant_property",
)

""" Build a deterministic, distinctive input set for the persistent-solver
    cross-call equivalence check. `seed` selects between two distinct
    starting topographies + sediment + parameter choices that exercise
    different dynamics inside the same toponum.
"""
function build_inputs(seed::Int)
    toponum = 1001
    xsize = 500_000.0
    topo_gridx = collect(range(0.0, xsize, length=toponum))
    topo_gridy_initial = if seed == 1
        [1000.0 * sin(2π * x / xsize) for x in topo_gridx]
    else
        [2000.0 * cos(3π * x / xsize) + 500.0 * sin(7π * x / xsize)
         for x in topo_gridx]
    end
    sediment_thickness_initial = fill(seed == 1 ? 200.0 : 50.0, toponum)
    params = SedimentTransportParameters(
        m2_yr_to_m2_s(0.25),
        m_yr_to_m_s(1.0),
        1e-4,
        m2_yr_to_m2_s(1e2),
        2000.0,
        years_to_seconds(5_000.0),
        years_to_seconds(seed == 1 ? 1.0e6 : 1.5e6),
        0.4,
        1/2000.0,
    )
    return (
        basic_grid_dims = (0.0, xsize),
        topo_gridx = topo_gridx,
        topo_gridy_initial = topo_gridy_initial,
        sediment_thickness_initial = sediment_thickness_initial,
        params = params,
        y_sealevel = seed == 1 ? 0.0 : 100.0,
        pelagic_rate = mm_yr_to_m_s(seed == 1 ? 0.0 : 0.5),
    )
end

function build_solver(inputs; use_optimized_solver::Bool=true)
    return SedimentTransportSolverManager.SedimentTransportSolver(
        inputs.basic_grid_dims,
        inputs.topo_gridx,
        inputs.topo_gridy_initial,
        inputs.sediment_thickness_initial,
        inputs.params,
        inputs.y_sealevel,
        inputs.pelagic_rate;
        _SOLVER_KWARGS...,
        use_optimized_solver = use_optimized_solver,
    )
end

function run_fresh(inputs; use_optimized_solver::Bool=true)
    solver = build_solver(inputs; use_optimized_solver=use_optimized_solver)
    SedimentTransportSolverManager.run_sediment_transport_time_steps!(
        solver, nothing)
    return solver
end

""" Capture every observable solver-state field that the production code
    reads after `run_sediment_transport_time_steps!` returns, plus the
    full `collections` history when `use_collections=true`. Snapshots are
    deep copies — the source solver can be mutated afterwards without
    invalidating them.
"""
function snapshot(solver)
    return (
        topo_gridy = copy(solver.topo_gridy),
        topo_gridy_initial = copy(solver.topo_gridy_initial),
        sediment_thickness_initial = copy(solver.sediment_thickness_initial),
        sediment_thickness_total_decompacted =
            solver.sediment_thickness_total_decompacted === nothing ?
                nothing : copy(solver.sediment_thickness_total_decompacted),
        sediment_thickness_initial_compacted =
            solver.sediment_thickness_initial_compacted === nothing ?
                nothing : copy(solver.sediment_thickness_initial_compacted),
        compaction_displacement_max =
            solver.compaction_displacement_max === nothing ?
                nothing : copy(solver.compaction_displacement_max),
        nsteps = solver.nsteps,
        transport_timestep = solver.transport_timestep,
        topo_collection =
            Dict(k => copy(v) for (k, v) in solver.collections.topo_collection),
        divides_collection =
            Dict(k => copy(v) for (k, v) in solver.collections.divides_collection),
        water_depth_collection =
            Dict(k => copy(v) for (k, v) in solver.collections.water_depth_collection),
        basement_collection =
            Dict(k => copy(v) for (k, v) in solver.collections.basement_collection),
    )
end

reset_solver!(solver, inputs) = SedimentTransportSolverManager.reset!(
    solver,
    inputs.basic_grid_dims,
    inputs.topo_gridx,
    inputs.topo_gridy_initial,
    inputs.sediment_thickness_initial,
    inputs.params,
    inputs.y_sealevel,
    inputs.pelagic_rate;
    _SOLVER_KWARGS...,
)

""" Cross-call equivalence regression for the persistent
    `SedimentTransportSolver` lifecycle.

    Validates the invariant that `reset!` carries no semantic state
    across calls beyond preallocated scratch buffers: a solver reused via
    `reset!` must produce bit-identical output to a solver freshly
    constructed for the same inputs.

    Three sub-checks per solver flavor (legacy + optimized):
      1. First run on a persistent solver matches a fresh solver run on
         the same inputs.
      2. After `reset!`-ing to a different input set, the persistent
         solver matches a fresh solver run on those new inputs.
      3. After a second `reset!` back to the original inputs, the
         persistent solver matches the original fresh-A baseline.

    The third check exercises the case the SedimentTransportTest cannot
    reach: that residual state from a previous run (collections, output
    fields, transport_timestep) is fully cleared by `reset!`.

    Each sub-check compares snapshots field-by-field with `==`; any
    floating-point drift fails it.

    Returns a NamedTuple of booleans suitable for `@test` assertions.
"""
function run_test(; use_optimized_solver::Bool=true)
    inputs_A = build_inputs(1)
    inputs_B = build_inputs(2)

    snap_A_fresh = snapshot(run_fresh(inputs_A;
        use_optimized_solver=use_optimized_solver))
    snap_B_fresh = snapshot(run_fresh(inputs_B;
        use_optimized_solver=use_optimized_solver))

    solver = build_solver(inputs_A; use_optimized_solver=use_optimized_solver)

    SedimentTransportSolverManager.run_sediment_transport_time_steps!(
        solver, nothing)
    snap_first = snapshot(solver)

    reset_solver!(solver, inputs_B)
    SedimentTransportSolverManager.run_sediment_transport_time_steps!(
        solver, nothing)
    snap_after_reset_to_B = snapshot(solver)

    reset_solver!(solver, inputs_A)
    SedimentTransportSolverManager.run_sediment_transport_time_steps!(
        solver, nothing)
    snap_after_reset_back_to_A = snapshot(solver)

    return (
        first_run_matches_fresh_A         = snap_first             == snap_A_fresh,
        after_reset_matches_fresh_B       = snap_after_reset_to_B  == snap_B_fresh,
        after_second_reset_matches_fresh_A = snap_after_reset_back_to_A == snap_A_fresh,
    )
end

""" Unit test for `SedimentTransportSolverManager.is_compatible`.

    `is_compatible` is the cache-validity gate inside `get_or_init_solver!`:
    a `false` return forces a full reconstruction of the persistent
    solver. Wrong logic here would either silently reuse a buffer of the
    wrong size (corruption / OOB) or always rebuild (defeating the
    hoist). All three branches need coverage:

    1. Matching `toponum` and `use_optimized_solver` → `true`.
    2. Mismatched `toponum` (any non-equal value) → `false`.
    3. Mismatched `use_optimized_solver` (the legacy/optimized flavors
       allocate different buffer sets, so they aren't interchangeable
       without reallocation) → `false`.

    Returns a NamedTuple of booleans that should all be `true` when
    `is_compatible` is correctly implemented.
"""
function run_is_compatible_test(; use_optimized_solver::Bool=true)
    inputs = build_inputs(1)
    solver = build_solver(inputs; use_optimized_solver=use_optimized_solver)
    toponum = length(inputs.topo_gridx)

    return (
        accepts_matching_inputs =
            SedimentTransportSolverManager.is_compatible(
                solver, toponum, use_optimized_solver) == true,
        rejects_toponum_mismatch =
            SedimentTransportSolverManager.is_compatible(
                solver, toponum + 1, use_optimized_solver) == false,
        rejects_optimized_flag_mismatch =
            SedimentTransportSolverManager.is_compatible(
                solver, toponum, !use_optimized_solver) == false,
    )
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    SedimentTransportSolverResetTest.run_test()
end
