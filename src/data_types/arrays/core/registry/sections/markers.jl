function get_markers_arrays()::NamedTuple
    return @arrays (
        #*************
        # Markers Arrays
        #*************

        # Grid Marker Relationship Arrays
        marker_xn = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Horizontal index of upper left node of basic grid cell containing marker.",
        ),
        marker_yn = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Vertical index of upper left node of basic grid cell containing marker.",
        ),
        marker_dx = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Normalized x-distance between marker and upper-left basic grid node in meters.",
        ),
        marker_dy = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Normalized y-distance between marker and upper-left basic grid node in meters.",
        ),
        marker_xn_vy = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Horizontal index of upper left node of staggered vy grid cell containing marker.",
        ),
        marker_yn_vy = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Vertical index of upper left node of staggered vy grid cell containing marker.",
        ),
        marker_dx_vy = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Normalized x-distance between marker and upper-left grid node of staggered vy grid in meters.",
        ),
        marker_dy_vy = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Normalized y-distance between marker and upper-left grid node of staggered vy grid in meters.",
        ),

        # Location Arrays
        marker_x = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : x-coordinate of marker in meters.",
        ),
        marker_y = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : y-coordinate of marker in meters.",
        ),

        # Material Arrays
        marker_matid = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Material ID of marker.",
        ),
        marker_rho = ArrayData(
            "kg/m^3", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker density in kg/m^3.",
        ),
        marker_porosity_initial = ArrayData(
            "fraction", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Initial porosity at sediment-water interface used in Athy's law in fraction.",
        ),
        marker_decay_depth = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Porosity decay depth used in Athy's law in meters.",
        ),
        marker_max_burial_depth = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Maximum burial depth of marker in meters.",
        ),
        marker_serpentinization = ArrayData(
            "fraction", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Serpentinization fraction of marker in fraction.",
        ),
        marker_serpentinization_heat_production = ArrayData(
            "W/m^3", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Serpentinization heat production of marker in W/m^3.",
        ),

        # Melt Arrays
        marker_meltfrac = ArrayData(
            "fraction", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Melt fraction of marker in fraction.",
        ),
        marker_extracted_meltfrac = ArrayData(
            "fraction", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Extracted melt fraction of marker in fraction.",
        ),
        marker_extractable_meltfrac = ArrayData(
            "fraction", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Extractable melt fraction of marker in fraction.",
        ),

        # Pressure Arrays
        marker_pr = ArrayData(
            "Pa", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pressure of marker in Pa.",
        ),

        # Rheology Arrays
        marker_eta = ArrayData(
            "Pa.s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker viscoplastic viscosity in Pa.s.",
        ),
        marker_fric_ini = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Initial friction coefficient (sine of friction angle) of marker.",
        ),
        marker_fric = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Friction coefficient (sine of friction angle) of marker.",
        ),
        marker_cohesion = ArrayData(
            "Pa", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Cohesion of marker in Pa.",
        ),
        marker_preexp = ArrayData(
            "1/s/MPa^n", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-exponential factor for dislocation creep of marker in 1/s/MPa^n.",
        ),
        marker_pfailure = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Flag indicating plastic failure of marker.",
        ),
        marker_eta_flow = ArrayData(
            "Pa.s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker viscosity from flow law in Pa.s.",
        ),
        marker_dilatation_angle = ArrayData(
            "Degrees", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker dilatation angle in Degrees.",
        ),

        # Strain Arrays
        marker_GII = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Accumulated strain of marker.",
        ),
        marker_strain_plastic = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Accumulated plastic strain of marker.",
        ),
        marker_exx = ArrayData(
            "1/s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Normal strain rate of marker in 1/s.",
        ),
        marker_exy = ArrayData(
            "1/s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Shear strain rate of marker in 1/s.",
        ),
        marker_sr_ratio = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Ratio of nodal-stress-based to nodal-velocity-based strain rate invariant for each marker.",
        ),
        marker_strain_rate_plastic = ArrayData(
            "1/s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Plastic strain rate of marker in 1/s.",
        ),
        marker_melt_damage = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Damage factor of marker associated with melt intrusion.",
        ),

        # Stratigraphy Arrays
        marker_age = ArrayData(
            "Myr", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Age of marker in million years.",
        ),

        # Stress Arrays
        marker_sxx = ArrayData(
            "Pa", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Deviatoric normal stress of marker in Pa.",
        ),
        marker_sxy = ArrayData(
            "Pa", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Deviatoric shear stress of marker in Pa.",
        ),

        # Thermal Arrays
        marker_TK = ArrayData(
            "K", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Temperature of marker in Kelvin.",
        ),
        marker_rhocp = ArrayData(
            "J/K/m^3", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker density multiplied by heat capacity in J/K/m^3.",
        ),
        marker_kt = ArrayData(
            "W/m/K", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker thermal conductivity in W/m/K.",
        ),
        marker_ha = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker adiabatic heating term (expansivity x temperature).",
        ),

        # Advection velocity/spin interpolated to markers (pre-allocated,
        # re-written every timestep by the Runge-Kutta interpolator).
        marker_vx = ArrayData(
            "m/s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Runge-Kutta-interpolated x-velocity at marker (m/s). "
            * "Re-computed every advection step; zero for markers outside the domain.",
        ),
        marker_vy = ArrayData(
            "m/s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Runge-Kutta-interpolated y-velocity at marker (m/s). "
            * "Re-computed every advection step; zero for markers outside the domain.",
        ),
        marker_spin = ArrayData(
            "1/s", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Runge-Kutta-interpolated spin rate at marker (1/s). "
            * "Re-computed every advection step; zero for markers outside the domain.",
        ),

        # Compaction scratch buffers used internally by
        # MarkerCompaction.compact_sediment_and_advect_markers and its
        # helpers (sticky path included). Each is single-use within one
        # compaction call (no cross-callsite sharing). Production callers
        # extract from model.markers.arrays.compaction.X.array and pass
        # explicitly to the function; the function fills them as needed.
        markers_topo_xindex_buffer = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for the topography-grid "
            * "x-index assigned to each sedimentary-basin marker by "
            * "MarkerCompaction.assign_topo_xindex_to_markers.",
        ),
        markers_compaction_yindex_buffer = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for the compaction-grid "
            * "y-index assigned to each sedimentary-basin marker by "
            * "MarkerCompaction.assign_compaction_yindex_to_markers.",
        ),
        markers_unit_distance_from_cell_top_buffer = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for the unit distance from "
            * "compaction-cell top per sedimentary-basin marker, "
            * "co-populated with markers_compaction_yindex_buffer.",
        ),
        total_marker_compaction_displacement_buffer = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for cumulative compaction "
            * "displacement per sedimentary-basin marker, populated by "
            * "MarkerCompaction.calculate_total_marker_compaction_displacement.",
        ),
        sticky_displacement_factors_buffer = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for displacement factors "
            * "computed by "
            * "MarkerCompaction.calculate_sticky_compaction_displacement_factors_opt.",
        ),
        sticky_marker_displacement_buffer = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for sticky-marker "
            * "displacements computed by "
            * "MarkerCompaction.calculate_sticky_marker_displacement_opt.",
        ),

        # Swarm-index scratch buffers used by
        # MarkerCompaction.calculate_swarm_indices_for_sediment_and_sticky and
        # MarkerCompaction.calculate_x_sorted_swarm_indices. Each is single-use
        # within one call and fully overwritten before the result is read, so
        # safe to reuse across the back-to-back sedimentary-basin and sticky
        # passes. The packed-prefix idiom matches
        # GridFuncs.get_indices_of_markers_outside_domain: fill in marker
        # order, then return a fresh `scratch[1:n_swarm]` so downstream
        # consumers continue to receive a tight `Vector{Int64}`.
        marker_swarm_index_scratch = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for the marknum-sized "
            * "intermediate inside "
            * "MarkerCompaction.calculate_swarm_indices_for_sediment_and_sticky. "
            * "Filled in marker-index order with matched marker indices "
            * "packed in front, then sliced into a tight return vector.",
        ),
        marker_swarm_x_gather_scratch = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for the gathered "
            * "`marker_x[swarm_indices]` values inside "
            * "MarkerCompaction.calculate_x_sorted_swarm_indices. Refilled "
            * "from `marker_x` per call before the sortperm! pass.",
        ),
        marker_swarm_sortperm_scratch = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Pre-allocated scratch holding the `sortperm!` "
            * "permutation inside "
            * "MarkerCompaction.calculate_x_sorted_swarm_indices. Initialized "
            * "by sortperm! itself per call (no carry-over).",
        ),

        # Solidification / shared random scratch used by Solidification.solidify!
        # and (via shared reuse) MarkerRecycle.RandomMarkerArray.get_random_marker_array.
        # Refilled via Random.rand! at each consumer call site.
        marker_random_buffer = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for per-marker random "
            * "numbers. Used by Solidification.solidify! during friction-"
            * "coefficient randomization, and by "
            * "MarkerRecycle.RandomMarkerArray.get_random_marker_array for "
            * "subsurface-marker reset. Each consumer refills via Random.rand! "
            * "before reading, so values do not leak between consumers.",
        ),

        # Marker recycling scratch buffer used by
        # GridFuncs.get_indices_of_markers_outside_domain to mark and pack
        # outside-marker indices in marker-index order.
        marker_outside_indices_scratch = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for the marknum-sized "
            * "intermediate inside "
            * "GridFuncs.get_indices_of_markers_outside_domain. Filled in "
            * "marker-index order, then packed in place before the function "
            * "returns a length-`nrecycle` copy.",
        ),
        marker_inside_flags_buffer = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Pre-allocated `Vector{Int8}` of length `marknum` "
            * "filled by `GridFuncs.get_marker_inside_flags` with `1` for "
            * "in-domain markers and `-1` for out-of-domain markers. The "
            * "function returns this buffer directly; downstream consumers "
            * "in the time loop treat it as read-only. Two call sites per "
            * "timestep (pre-solver setup and post-solver re-evaluation) "
            * "overwrite the buffer in sequence — the second call's values "
            * "supersede the first; nothing in between persists a reference "
            * "to the prior values past the second call.",
        ),

        # Serpentinization scratch buffer used by
        # Serpentinization.calculate_marker_serpentinization.
        marker_serpentinization_increment_buffer = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for the per-marker "
            * "incremental serpentinization ratio. Refilled per call to "
            * "Serpentinization.calculate_marker_serpentinization.",
        ),

        # Subgrid heat diffusion scratch buffer used by
        # SubgridDiffusion.calculate_subgrid_temperature_change_and_correct_marker_temperature!.
        marker_subgrid_temp_delta_buffer = ArrayData(
            "K", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for the per-marker subgrid "
            * "temperature change. Refilled per call to "
            * "SubgridDiffusion.calculate_subgrid_temperature_change_and_correct_marker_temperature!.",
        ),

        # Lithostatic pressure column-filter scratch buffers used by
        # LithostaticPressure.filter_markers_for_column. Three marknum-sized
        # buffers each pack column-matching markers in front; the function
        # then copies the prefix into 3 small tight return vectors. Each
        # buffer is single-use within filter_markers_for_column and never
        # shared with any other call site or cross-domain consumer.
        marker_x_filter_scratch = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated tmp buffer for marker x-coordinates "
            * "during column filtering inside "
            * "LithostaticPressure.filter_markers_for_column.",
        ),
        marker_y_filter_scratch = ArrayData(
            "m", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated tmp buffer for marker y-coordinates "
            * "during column filtering inside "
            * "LithostaticPressure.filter_markers_for_column.",
        ),
        marker_rho_filter_scratch = ArrayData(
            "kg/m^3", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated tmp buffer for marker densities "
            * "during column filtering inside "
            * "LithostaticPressure.filter_markers_for_column.",
        ),

        # Sediment-transport marker advection scratch buffers used by the four
        # functions in MarkerAdvection.jl
        # (calculate_sediment_compaction_displacement_factors!,
        # calculate_sediment_marker_displacement!,
        # calculate_sticky_compaction_displacement_factors!,
        # calculate_sticky_marker_displacement!). The two buffers are reused
        # across sediment and sticky advection paths within
        # advect_markers_using_compaction; each writer fill!(0.0)s before
        # writing, so prior values never leak between consumers.
        sediment_transport_marker_displacement_factors_buffer = ArrayData(
            "None",
            MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for displacement factors "
            * "computed by MarkerAdvection.calculate_sediment_compaction_displacement_factors! "
            * "and MarkerAdvection.calculate_sticky_compaction_displacement_factors!. "
            * "Filled with fill!(0.0) at start of each writer call.",
        ),
        sediment_transport_marker_displacement_buffer = ArrayData(
            "m",
            MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Pre-allocated scratch for marker displacements "
            * "computed by MarkerAdvection.calculate_sediment_marker_displacement! "
            * "and MarkerAdvection.calculate_sticky_marker_displacement!. "
            * "Filled with fill!(0.0) at start of each writer call.",
        ),
    )
end
