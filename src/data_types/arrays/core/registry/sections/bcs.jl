function get_bcs_arrays()::NamedTuple
    return @arrays (
        #*************
        # Boundary Conditions Arrays
        #*************

        # Temperature Boundary Conditions Arrays
        btopt = ArrayData(
            "K", BcArrayFloatState, "one-dimensional",
            "`(xnum, 2)` : Temperature boundary condition along top boundary: \n"
            *"\t tk[1][j] = btopt[j][1] + tk[2][j]*btop[j][2]",
        ),
        bbottomt = ArrayData(
            "K", BcArrayFloatState, "one-dimensional",
            "`(xnum, 2)` : Temperature boundary condition along bottom boundary: \n"
            *"\t tk[ynum][j] = bbottomt[j][1] + tk[ynum-1][j]*bbottomt[j][2]",
        ),
        bleftt = ArrayData(
            "K", BcArrayFloatState, "one-dimensional",
            "`(ynum, 2)` : Temperature boundary condition along left boundary: \n"
            *"\t tk[i][1] = bleftt[i][1] + bleftt[i][2]*tk[i][2]",
        ),
        brightt = ArrayData(
            "K", BcArrayFloatState, "one-dimensional",
            "`(ynum, 2)` : Temperature boundary condition along right boundary: \n"
            *"\t tk[i][xnum] = brightt[i][1] + brightt[i][2]*tk[i][xnum-1]",
        ),

        # Component-wise Velocity Boundary Conditions Arrays
        btopx = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vx boundary condition along top boundary: \n"
            *"\t vx[1,j,:] = btopx[j,1] + vx[2,j,:]*btopx[j,2]",
        ),
        btopy = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vy boundary condition along top boundary: \n"
            *"\t vy[1,j,:] = btopy[j,1] + vy[2,j,:]*btopy[j,2]",
        ),
        btopz = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vz boundary condition along top boundary: \n"
            *"\t vz[1,j,:] = btopz[j,1] + vz[2,j,:]*btopz[j,2]",
        ),
        bbottomx = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vx boundary condition along bottom boundary: \n"
            *"\t vx[ynum+1,j,:] = bbottomx[j,1] + vx[ynum, j,:]*bbottomx[j,2]",
        ),
        bbottomy = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vy along bottom boundary: \n"
            *"\t vy[ynum,j,:] = bbottomy[j,1] + vy[ynum-1, j,:]*bbottomy[j,2]",
        ),
        bbottomz = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 2)` : Vz boundary condition along bottom boundary: \n"
            *"\t vz[ynum+1,j,:] = bbottomz[j,1] + vz[ynum, j,:]*bbottomz[j,2]",
        ),
        bleftx = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vx boundary condition along left boundary: \n"
            *"\t vx[i,1] = bleftx[i,1] + vx[i, 2]*bleftx[i,2]",
        ),
        blefty = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vy boundary condition along left boundary: \n"
            *"\t vy[i,1] = blefty[i,1] + vy[i, 2]*blefty[i,2]",
        ),
        bleftz = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vz boundary condition along left boundary: \n"
            *"\t vz[i,1,:] = bleftz[i,1] + vz[i,2,:]*bleftz[i,2]",
        ),
        brightx = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vx boundary condition along right boundary: \n"
            *"\t vx[i,xnum] = brightx[i, 1] + vx[i, xnum-1]*brightx[i,2]",
        ),
        brighty = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vy boundary condition along right boundary: \n"
            *"\t vy[i,xnum+1] = brighty[i,1] + vx[i, xnum]*brighty[i,2]",
        ),
        brightz = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vz boundary condition along right boundary: \n"
            *"\t vz[i,xnum+1,:] = brightz[i,1] + vz[i,xnum,:]*brightz[i,2]",
        ),
        bfrontx = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vx boundary condition along front boundary: \n"
            *"\t vx[i,j,1] = bfrontx[i,j,1] + vx[i,j,2]*bfrontx[i,j,2]",
        ),
        bfronty = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vy boundary condition along front boundary: \n"
            *"\t vy[i,j,1] = bfronty[i,j,1] + vy[i,j,2]*bfronty[i,j,2]",
        ),
        bfrontz = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vz boundary condition along front boundary: \n"
            *"\t vz[i,j,1] = bfrontz[i,j,1] + vz[i,j,2]*bfrontz[i,j,2]",
        ),
        bbackx = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vx boundary condition along back boundary: \n"
            *"\t vx[i,j,znum+1] = bbackx[i,j,1] + vx[i,j,znum]*bbackx[i,j,2]",
        ),
        bbacky = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vy boundary condition along back boundary: \n"
            *"\t vy[i,j,znum+1] = bbacky[i,j,1] + vy[i,j,znum]*bbacky[i,j,2]",
        ),
        bbackz = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 2)` : Vz boundary condition along back boundary: \n"
            *"\t vz[i,j,znum  ] = bbackz[i,j,1] + vz[i,j,znum-1]*bbackz[i,j,2]",
        ),

        # Combined Velocity Boundary Conditions Arrays
        btop = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 6)` : Velocity boundary conditions along top boundary: \n"
            *"\t vx[1,j,:] = btop[j,1] + vx[2,j,:]*btop[j,2] \n"
            *"\t vy[1,j,:] = btop[j,3] + vy[2,j,:]*btop[j,4] \n"
            *"\t vz[1,j,:] = btop[j,5] + vz[2,j,:]*btop[j,6]"
        ),
        bbottom = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(xnum+1, 6)` : Velocity boundary conditions along bottom boundary: \n"
            *"\t vx[ynum+1,j,:] = bbottom[j,1] + vx[ynum  ,j,:]*bbottom[j,2] \n"
            *"\t vy[ynum  ,j,:] = bbottom[j,3] + vy[ynum-1,j,:]*bbottom[j,4] \n"
            *"\t vz[ynum+1,j,:] = bbottom[j,5] + vz[ynum  ,j,:]*bbottom[j,6]",
        ),
        bleft = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 6)` : Velocity boundary conditions along left boundary: \n"
            *"\t vx[i,1,:] = bleft[i,1] + vx[i,2,:]*bleft[i,2] \n"
            *"\t vy[i,1,:] = bleft[i,3] + vy[i,2,:]*bleft[i,4] \n"
            *"\t vz[i,1,:] = bleft[i,5] + vz[i,2,:]*bleft[i,6]",
        ),
        bright = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 6)` : Velocity boundary conditions along right boundary: \n"
            *"\t vx[i,xnum  ,:] = bright[i,1] + vx[i,xnum-1,:]*bright[i,2] \n"
            *"\t vy[i,xnum+1,:] = bright[i,3] + vy[i,xnum  ,:]*bright[i,4] \n"
            *"\t vz[i,xnum+1,:] = bright[i,5] + vz[i,xnum  ,:]*bright[i,6]",
        ),
        bfront = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 6)` : Velocity boundary conditions along front boundary: \n"
            *"\t vx[i,j,1] = bfront[i,1] + vx[i,j,2]*bfront[i,2] \n"
            *"\t vy[i,j,1] = bfront[i,3] + vy[i,j,2]*bfront[i,4] \n"
            *"\t vz[i,j,1] = bfront[i,5] + vz[i,j,2]*bfront[i,6]",
        ),
        bback = ArrayData(
            "m/s", BcArrayFloatState, "one-dimensional",
            "`(ynum+1, 6)` : Velocity boundary conditions along back boundary: \n"
            *"\t vx[i,j,znum+1] = bback[i,1] + vx[i,j,znum]*bback[i,2] \n"
            *"\t vy[i,j,znum+1] = bback[i,3] + vy[i,j,znum]*bback[i,4] \n"
            *"\t vz[i,j,znum  ] = bback[i,5] + vz[i,j,znum-1]*bback[i,6]",
        ),
    )
end
