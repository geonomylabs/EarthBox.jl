function get_interpolation_arrays()::NamedTuple
    return @arrays (
        #*************
        # Interpolation Arrays
        #*************

        # Grid Weights Arrays
        wtnodes = ArrayData(
            "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Summed marker weight factors on basic grid nodes "
            * "using inclusive method with maximum search radius for markers.",
        ),
        wtetas = ArrayData(
            "None", ScalarArray2DState, "basic",
            "`(ynum, xnum)` : Summed marker weight factors on basic grid nodes "
            * "using exclusive method with shorter search radius for markers.",
        ),
        wtetan = ArrayData(
            "None", ScalarArray2DState, "pressure",
            "`(ynum - 1, xnum - 1)` : Summed marker weight factors for normal "
            * "viscosity on pressure nodes.",
        ),
        wtnodes_vy = ArrayData(
            "None", ScalarArray2DState, "vy",
            "`(ynum + 1, xnum)` : Summed marker weight factors on Vy staggered "
            * "grid nodes using inclusive method with maximum search radius.",
        ),

        # Marker Weights Arrays
        marker_wtforULnode = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for upper-left basic grid node.",
        ),
        marker_wtforLLnode = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for lower-left basic grid node.",
        ),
        marker_wtforURnode = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for upper-right basic grid node.",
        ),
        marker_wtforLRnode = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for lower-right basic grid node.",
        ),
        marker_wtforCnode = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for central pressure grid node.",
        ),
        marker_wtforULnodeVy = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for upper-left Vy staggered grid node.",
        ),
        marker_wtforLLnodeVy = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for lower-left Vy staggered grid node.",
        ),
        marker_wtforURnodeVy = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for upper-right Vy staggered grid node.",
        ),
        marker_wtforLRnodeVy = ArrayData(
            "None", MarkerArrayFloat1DState, "NA",
            "`(marknum)` : Marker weight for lower-right Vy staggered grid node.",
        ),
    )
end
