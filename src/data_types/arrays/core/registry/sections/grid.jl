function get_grid_arrays()::NamedTuple
    return @arrays (
        #*************
        # Grids Arrays
        #*************

        # Basic Grid Arrays
        gridx_b = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(xnum)` : x-coordinates of basic grid nodes in meters.",
        ),
        xstp_b = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(xnum - 1)` : Width of cells in x-direction for basic grid in meters.",
        ),
        gridy_b = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(ynum)` : y-coordinates of basic grid nodes in meters.",
        ),
        ystp_b = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(ynum - 1)` : Width of cells in y-direction for basic grid in meters.",
        ),
        # Pressure Grid Arrays
        gridy_pr = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(ynum - 1)` : y-locations of pressure grid nodes in meters.",
        ),
        ystp_pr = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(ynum - 2)` : Width of pressure cells in y-direction in meters.",
        ),
        gridx_pr = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(xnum - 1)` : x-locations of pressure grid nodes in meters.",
        ),
        xstp_pr = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(xnum - 2)` : Width of pressure cells in x-direction in meters.",
        ),
        # Staggered Vx Grid Arrays
        gridy_vx = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(ynum + 1)` : y-locations of staggered Vx grid nodes including ghost nodes in meters.",
        ),
        ystp_vx = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(ynum)` : Width of staggered vx grid cells in y-direction in meters.",
        ),
        # Staggered Vy Grid Arrays
        gridx_vy = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(xnum + 1)` : x-locations of staggered Vy grid nodes including ghost nodes in meters.",
        ),
        xstp_vy = ArrayData(
            "m", GridArray1DState, "one-dimensional",
            "`(xnum)` : Width of staggered vy grid cells in x-direction in meters.",
        ),
    )

end
