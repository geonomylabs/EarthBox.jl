function get_melting_arrays()::NamedTuple
    return @arrays (
        #*************
        # Melting Arrays
        #*************

        # Extraction Arrays
        xstart_drainage = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Current x-locations in meters of the start of the drainage area. "
            * "The drainage area is the region where the melt migrates to the "
            * "shallowest position.",
        ),
        xend_drainage = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Current x-locations in meters of the end of the drainage area. "
            * "The drainage area is the region where the melt migrates to the "
            * "shallowest position.",
        ),
        melt_residuals = ArrayData(
            "m^3", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Residual melt volume in m^3 remaining in the extrusion volume "
            * "after the melt has been extracted in drainage basin.",
        ),
        extrusion_volumes = ArrayData(
            "m^3", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Volume of extracted melt in m^3 that will be extruded to surface.",
        ),
        xmid_molten_zones = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : X-locations in meters of the midpoints of the molten region.",
        ),
        ytop_molten_zones = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Y-locations in meters of the top of the molten region.",
        ),
        width_molten_zones = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Width of the molten region in drainage basin in meters.",
        ),
        avg_shallow_partial_melt_xcoors = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Average x-coordinates of the shallowest partial melt "
            * "marker in drainage basins in meters.",
        ),
        avg_shallow_partial_melt_ycoors = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Average y-coordinates of the shallowest partial melt "
            * "marker in drainage basins in meters.",
        ),
        xstart_drainage_o = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Old x-locations of the start of the drainage area. "
            * "The drainage area is the region where the melt migrates to the "
            * "shallowest position in meters.",
        ),
        xend_drainage_o = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(nmax_basin)` : Old x-locations of the end of the drainage area. "
            * "The drainage area is the region where the melt migrates to the "
            * "shallowest position in meters.",
        ),

        # Extraction scratch buffers (marker-sized, reused across timesteps)
        partial_melt_marker_indices = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Pre-allocated scratch buffer for packed indices of "
            * "partially molten markers. Only positions [1:nmarkers_partial_melt] "
            * "are valid after each call; the tail is stale. Initialized to 0.",
        ),
        marker_indices_tmp = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Pre-allocated scratch buffer for packed indices of "
            * "markers in the mantle injection search domain. Only positions "
            * "[1:nmarkers_injection_domain] are valid after each call; the tail "
            * "is stale. Initialized to -1.",
        ),
        pm_layer_counts = ArrayData(
            "None", Array1DIntState, "NA",
            "`(nlayers)` : Pre-allocated scratch buffer for per-layer counts of "
            * "partially molten mantle markers used by the layered "
            * "shallowest-marker search. Populated by "
            * "construct_layered_partially_molten_arrays!.",
        ),
        pm_layer_offsets = ArrayData(
            "None", Array1DIntState, "NA",
            "`(nlayers + 1)` : Pre-allocated scratch buffer holding exclusive "
            * "prefix-sum offsets into layered_partial_melt_indices. Layer i's "
            * "valid indices occupy positions [offsets[i]+1 : offsets[i+1]].",
        ),
        layered_partial_melt_indices = ArrayData(
            "None", MarkerArrayInt1DState, "NA",
            "`(marknum)` : Pre-allocated scratch buffer holding partially "
            * "molten mantle marker indices packed by layer (counting-sort "
            * "style). Use pm_layer_counts and pm_layer_offsets to bound per-"
            * "layer reads.",
        ),
    )
end
