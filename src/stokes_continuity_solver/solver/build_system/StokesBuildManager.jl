""" Build compact system of equations for Stokes-continuity equations.

Description of the 2D Staggered Grid
================================================================================

::

        j   1            2            3             4           5            6

          xstp_vy       xstpc2       xstpc3       xstpc4       xstpc5       xstpc6
     |............|............|............|............|............|............|
                xstp_b        xstp2        xstp3        xstp4        xstp5        xstp6
            |............|............|............|............|............|............|

            vx           vx           vx           vx           vx           vx             .....
    i

    1 vy    +-----vy-----+-----vy-----+-----vy-----+-----vy-----+-----vy-----+     vy ..... ytspc1
            |1           |5           |9           |13          |17          |
            |            |            |            |            |            |
            vx   p'3    vx1   p'15   vx13   p'27  vx25  p'39   vx37  p'51   vx49      ystp0 .....
            |            |            |            |            |            |
            |            |            |            |            |            |
    2 vy    +-----vy2----+----vy14----+----vy26----+----vy38----+----vy50----+     vy ..... ytspc2
            |2           |6           |10          |14          |18          |
            |            |            |            |            |            |
            vx    p'6   vx4   p'18   vx16   p'30  vx28  p'42   vx40  p'54   vx52      ystp_b .....
            |            |            |            |            |            |
            |            |            |            |            |            |
    3 vy    +-----vy5----+----vy17----+----vy29----+----vy41----+----vy53----+     vy ..... ytspc3
            |3           |7           |11          |15          |19          |
            |            |            |            |            |            |
            vx   p'9    vx7   p'21   vx19   p'33  vx31  p'45   vx43  p'57   vx55      ystp2 .....
            |            |            |            |            |            |
            |            |            |            |            |            |
    4 vy     +-----vy8----+----vy20----+----vy32----+----vy44----+----vy56----+     vy ..... ytspc4
            |4           |8           |12          |16          |20          |
            |            |            |            |            |            |
            vx   p'12   vx10  p'24   vx22   p'36  vx34  p'48   vx46  p'60   vx58      ystp3 .....
            |            |            |            |            |            |
            |            |            |            |            |            |
    5 vy     +----vy11----+----vy23----+----vy35----+----vy47----+----vy59----+     vy ..... ytspc5


            vx           vx           vx           vx           vx           vx             .....

Figure 1: 2D staggered grid and numbering scheme used for discretizing
Stokes and continuity equations. Indices i and j apply to basic grid
nodes denoted with "+". Symbols vx and vy refer to staggered sub-grid
nodes for the x- and y-components of velocity respectively. p' refers
to the staggered sub-grid nodes for scaled pressure with nodes located
in the center of basic grid cells. Note that staggered sub-grids have
different indices than the ones displayed in this figure for the basic
grid.

The staggered grid is composed of a basic grid and 3 staggered sub-grids
(Figure 1). The basic grid has dimensions (xnum, ynum). Basic grid nodes
are denoted with "+" in Figure 1. Shear stress (sxy), shear viscosity
(etas), shear modulus (mus) and density (rho) are defined at basic nodes.
When discretizing the heat equation temperature (tk), thermal conductivity
(kt) and radiogenic heat production are also defined on the basic grid.

The staggered sub-grids include (1) a grid with dimensions (xnum, ynum+1)
where the x-component of velocity vx is defined, (2) a grid with
dimensions (xnum+1, ynum) where the y-component of velocity vy is defined
and (3) a grid with dimensions (xnum-1, ynum-1) where scaled pressure p'
and normal viscosity (etan) are defined. vx and vy nodes located along
boundaries and outside of the basic grid are referred to as the ghost
boundary nodes that are used to define boundary conditions. The ghost
boundary nodes with numbers are ghost unknowns used as place holders in
the sparse matrix of the system of equations described below.

The x-spacing and y-spacing of the basic grid cells are stored in arrays
xstp(xnum-1) and ystp(ynum-1). The x-spacing of the staggered vy grid
centered on basic nodes is stored in xstpc(xnum). The y-spacing of the
staggered vx grid centered on basic nodes is stored in ystpc(ynum). See
Figure 1 for details.

The following arrays are used to store information defined on the basic
grid:

    rho(ynum, xnum)  : density (kg/m^3)
    sxy(ynum, xnum)  : shear stress (Pa)
    etas(ynum, xnum) : viscosity for shear stress (Pa.s)
    mus(ynum, xnum)  : shear modulus (Pa)
    tk(ynum, xnum)   : temperature (K)
    kt(ynum, xnum)   : thermal conductivity (W/m/K)
    hr(ynum, xnum)   : Radiogenic heat production (W/m^3)

The following arrays are used to store information on the staggered grids
for the x-component and y-component of velocity:

    Vx(ynum+1, xnum) : x-component of velocity (m/s)
    Vy(ynum, xnum+1) : y-component of velocity (m/s)

The following arrays are used to store information on the staggered grid
for pressure:

    p(xnum-1, ynum-1)    : pressure (Pa)
    p'(xnum-1, ynum-1)   : scaled pressure
    sxx(xnum-1, ynum-1)  : normal stress (Pa)
    etan(xnum-1, ynum-1) : viscosity for normal stress (Pa.s)

Unscaled pressure p is related to scaled pressure p' according to the
following equation:

    p = pscale*p'                                                    (Eq. 1)

where pscale is equal to

    pscale = 2.0*etan[0,0]/(xstpavr + ystpavr)                       (Eq. 2)

and etan[0,0] is the normal viscosity at the first basic node, xstpavr is
the average grid spacing in the x-direction, and ystpavr is the average
the average spacing in the y-direction.


Discretizing Governing Equations Using Finite Difference Stencils
================================================================================

Coefficients and right-hand-side terms are generated by looping over
solution nodes of the basic grid (i.e. grid cells) and discretizing the
governing equations associated with unknowns vx, vy and p' using the
relevant stencil.

::

        x-Stokes Stencil applied for each vx unknown
        +------------------------+------------------------+
        |                        |                        |
        |                        |                        |
        |                        |                        |
        |                        |                        |
        |                     Vx(i,j+1)                   |---------
        |                       vxU                       |
        |                        |                        |
        |                        |                        |
        |                        |                        |
     BN(i,j)------Vy(i,j+1)---etas(i,j+1)---Vy(i,j+2)-----+  ystpc(i) -------
        |           vyUL        etasU        vyUR         |    dyU
        |                       rhoU                      |
        |                        |                        |
        |        P'(i,j)         |         P'(i,j+1)      |
    Vx(i+1,j)      p'L      vx(i+1,j+1)      p'R    Vx(i+1,j+2)----- ystp(i)
        vL       etan(i,j)      vxC       etan(i,j+1)    vxR          dyC
        |         etanL          |          etanR         |
        |                        |                        |
        |                        |                        |
        +-----Vy(i+1,j+1)---etas(i+1,j+1)--Vy(i+1,j+2)----+ ystpc(i+1)-------
        |         vyLL         etasD          vyLR        |    dyD
        |                      rhoD                       |
        |                        |                        |
        |                        |                        |
        |                    Vx(i+2,j+1)                  |---------
        |                       vxD                       |
        |                        |                        |
        |                        |                        |
        |                        |                        |
        +------------------------+------------------------+
        |         xstp(j)        |        xstp(j+1)       |
                    dxL                      dxR
                    |       xstpc(j+1)       |
                                dxC

Figure 2: x-Stokes equation stencil. Indices i and j refer to the basic
grid. Indices for nodes of staggered sub-grids Vx, Vy and P' are
calculated using i and j from basic nodes. BN refers to the basic node
associated with this stencil. Abbreviated terms are shown below array
elements that denote the location of nodes relative to the central vx
node. These abbreviated terms use the following suffixes: "C" refers to
the central vx location, "U" refers to a location up above
central node, "D" refers to a location down below the central node,
"L" refers to a location to the left of the central node, "R"
refers to a location to the right of the central node, "UL" refers to
a location to the upper left of the central node, UR refers to a
location to the upper right of the central node, "LL" refers to a
location the lower left of the central node and "LR" refers to a
location to the lower right of the central node.

   The following x-Stokes equation, which describes force balance per unit
   volume, is formulated at the central node of the stencil associated with
   the unknown vxC:

       dSIGMAxx/dx+dSIGMAxy/dy-dP/dx=RX                                  (Eq.3)

   where SIGMAxx/dx is normal stress gradient in the x-direction, SIGMAxy/dy
   is shear stress gradient in the y-direction and dP/dx is the pressure
   gradient in the x-direction and RX is the right-hand-side term. Using the
   stencil from Figure (2), the discretized form of (Eq. 3) is as
   follows:

         [ 2*etanR*(vxR - vxC)/dxR - 2*etanL*(vxC - vxL)/dxL ]/dxC
       + { etasD*[ (vxD - vxC)/dyD + (vyLR - vyLL)/dxC ]
          - etasU*[ (vxC - vxU)/dyU + (vyUR - vyUL)/dxC ] }/dyC
       - pscale*(p'R - p'L)/dxC = -gravity_x*(rhoU+rhoD)/2                     (Eq. 4)

   Rearranging (Eq. 4) leads to the following equation with coefficients for
   unknowns:

    - (2*etanL/dxL/dxC + 2*etanR/dxR/dxC + etasU/dyU/dyC + etasD/dyD/dyC)*vxC
    + 2*etanL/dxL/dxC*vxL
    + 2*etanR/dxR/dxC*vxR
    + etasU/dyU/dyC*vxU
    + etasD/dyD/dyC*vxD
    + etasU/dxC/dyC*vyUL
    - etasD/dxC/dyC*vyLL
    - etasU/dxC/dyC*vyUR
    - etasD/dxC/dyC*vyLR
    - pscale/dxC*p'R
    + pscale/dxC*p'L
    = -gravity_x*(rhoU+rhoD)/2                                                 (Eq. 7)

    In cases where ghost boundary nodes are included in the stencil formulation
    some unknowns become constants requiring regrouping of terms and the
    reformulation of some coefficients and the right-hand-side term.

    The x-Stokes stencil is not applied along basic cells located along the
    right boundary since the velocity for the central node is prescribed. For
    these cases the matrix element (ivx, ivx) is equal to 1 and the
    right-hand-side is equal to the prescribed vx velocity (vx = zero for free
    slip and no slip). These nodes are referred to as ghost unknowns by Gerya (2010).

y-Stokes Stencil applied for each vy unknown
+--------------------BN(i,j)----Vy(i,j+1)---------+------------------------+-----
|                       |          vyU            |                        |
|                       |                         |                        |
|                       |                         |                        |
|                       |        P'(i,j)          |                        |
|                   Vx(i+1,j)      p'U        Vx(i+1,j+1)                  |ystp(i)-----
|                      vxUL      etan(i,j)       vxUR                      | dyU
|                       |         etanU           |                        |
|                       |                         |                        |
|                       |                         |                        |
+----Vy(i+1,j)---etas(i+1,j)--Vy(i+1,j+1)--etas(i+1,j+1)---Vy(i+1,j+2)-----+--ystpc(i+1)
|        vyL          etasL        vyC          etasR          vyR         |   dyC
|                     rhoL                      rhoR                       |
|                       |       p'(i+1,j)         |                        |
|                       |          p'D            |                        |
|                   Vx(i+2,j)   etan(i+1,j)   Vx(i+2,j+1)                  |ystp(i+1)---
|                      vxLL       etanD          vxLR                      | dyD
|                       |                         |                        |
|                       |                         |                        |
|                       |                         |                        |
+-----------------------+------Vy(i+2,j+1)--------+------------------------+-----
                                   vyD

           |         xstpc(j)       |        xstpc(j+1)       |
                       dxL                      dxR
                        |         xstp(j)         |
                                   dxC
        Figure 3: y-Stokes equation stencil

   The following y-Stokes equation, which describes force balance per unit
   volume, is formulated at the central node of the stencil associated with
   the unknown vyC:

       dSIGMAyy/dy + dSIGMAyx/dx - dP/dy = RY                            (Eq.8)

   where SIGMAyy/dx is normal stress gradient in the y-direction, SIGMAyx/dx
   is shear stress gradient in the x-direction and dP/dy is the pressure
   gradient in the y-direction and RY is the right-hand-side term. Using the
   stencil from Figure (3), the discretized form of (Eq. 8) is as
   follows:

       [ 2*etanD*(vyD - vyC)/dyD - 2*etanU*(vyC - vyU)/dyU) ]/dyC
     + { etasR*[ (vyR - vyC)/dxR + (vxLR - vxUR)/dyC) ]
        - etasL*[ (vyC - vyL)/dxL + (vxLL - vxUL)/dyC ] }/dxC
     - pscale*(p'D - p'U)/dyC = -gravity_y*(rhoR+rhoL)/2                (Eq. 9)

   Rearranging (Eq. 9) leads to the following equation with coefficients for
   unknowns:

    - (2*etanU/dyC/dyU + 2*etanD/dyC/dyD + etasL/dxC/dxL + etasR/dxC/dxR)*vyC
    + 2*etasL/dxC/dxL*vyL
    + 2*etasR/dxC/dxR*vyR
    + 2*etanU/dyC/dyU*vyU
    + 2*etanD/dyC/dyD*vyD
    + etasL/dxC/dyC*vxUL
    - etasL/dxC/dyC*vxLL
    - etasR/dxC/dyC*vxUR
    - etasR/dxC/dyC*vxLR
    - pscale/dyC*p'D
    + pscale/dyC*p'U
    = -gravity_y*(rhoR+rhoL)/2                                                (Eq. 10)


                Continuity stencil applied for each p' unknown

                    BN(i,j)------Vy(i,j+1)------+        ......
                        |          vyU          |
                        |                       |
                        |                       |
                        |                       |
                    Vx(i+1,j)     P'(i,j)  Vx(i+1,j+1)   ystp(i)
                       vxL         p'C         vxR         dy
                        |                       |
                        |                       |
                        |                       |
                        |                       |
                        +------Vy(i+1,j+1)------+       ......
                                   vyD

                        |........xstp(j).......|
                                  dx
        Figure 4: Continuity equation stencil.

   The following continuity equation is formulated at the central node of the
   stencil associated with the unknown p':

       dvx/dx + dvy/dy = RC                                             (Eq.11)

    where RC is the right-hand-side term and other parameters are defined
    above. For cases where RC = 0 and using the stencil from Figure (4), the
    discretized form of (Eq. 11) is as follows:

        pscale*[ (vxR - vxL)/dx + (vyD - vyU)/dy] = 0                  (Eq. 12)

    Rearranging (Eq. 12) leads to the following equation with coefficients for
    unknowns:

          pscale/dx*vxR
        - pscale/dx*vxL
        + pscale/dy*vyD
        - pscale/dy*vyU
        = 0

    Building the System of Equations
    ============================================================================

    A system of equations is built with the following form:

        LX = R                                                         (Eq. 13)

    where L is a sparse matrix with size NxN, X is a vector of unknowns and R
    is a vector containing the right-hand-side terms. Note that the matrix L is
    also referred to as the large discretized matrix in this section. N is
    related to the basic grid nodes according to the following equation:

        N = (xnum-1)*(ynum-1)*3                                        (Eq. 14)

    where xnum is the number of basic grid nodes in the horizontal x-direction
    (j-index) and ynum is the number of basic nodes in the vertical
    y-direction (i-index). Conceptually this equation computes the number of
    basic grid cells, which is then multiplied by the number of unknowns per
    cells (i.e. 3).

    The general structure of L for a case with xnum = ynum = 4 is shown in
    Figure 5 where i' and j' refer to large matrix indices for rows and column
    respectively. Stencils for x-Stokes, y-Stokes and continuity are applied
    for each respective row i' leading to coefficients for a central node from
    the stencil corresponding to the unknown associated with the row. The
    central nodes have coefficients denoted as "c" in Figure 5 and are located
    at the matrix location (i',j') where i' = j'.

    Values for L are filled by looping over all columns of the basic grid
    (index j) associated with a given row (index i).

    The appropriate stencil is applied for each matrix row associated
    with an unknown of a particular type, Coefficients are defined for all
    unknowns associated with the stencil and are stored in rows of the large
    discretized matrix L at column locations associated with the particular
    unknown.

    A vector R is also filled during the looping process with right-hand-side
    terms associated with the stencil applied to generate coefficients of a
    given row of the matrix.

    The central node for the x-Stokes equation is located within the large
    matrix at (ivx, ivx) where

        cell_index  = (j-1)*(ynum-1) + i                                (Eq. 15)

        ivx = cell_index*3 - 2                                          (Eq. 16)

    and i and j refer to indices of the basic grid. cell_index is the index for
    the basic grid cells associated with (i,j) (see Figure 1).

    Coefficients for the central node of the y-Stokes stencil are located at
    matrix location (ivx+1, ivx+1), and the coefficients for the central
    node of the continuity equation are located at matrix location
    (ivx+2, ivx+2).

    When a stencil is applied to fill elements in a large-matrix row,
    coefficients are generated for additional unknowns. These coefficients are
    stored at their respective locations within a given large matrix row
    defined with respect to the central node.

    For the x-Stokes equation coefficient matrix locations are defined as
    follows:

        vxC  : (ivx, ivx)
        vxL  : (ivx, ivx-hshift_to_vxR)
        vxR  : (ivx, ivx+hshift_to_vxR)
        vxU  : (ivx, ivx-3)
        vxD  : (ivx, ivx+3)
        vyUL : (ivx, ivx-2)
        vyUR : (ivx, ivx-2+hshift_to_vxR)
        vyLL : (ivx, ivx+1)
        vyLR : (ivx, ivx+1+hshift_to_vxR)
        p'L  : (ivx, ivx+2)
        p'R  : (ivx, ivx+2+hshift_to_vxR)

    where

    hshift_to_vxR = (ynum-1)*3                                         (Eq. 17)

    For the y-Stokes equation coefficient matrix locations are defined as
    follows:

        vyC  : (ivy, ivy)
        vyU  : (ivy, ivy-3)
        vyD  : (ivy, ivy+3)
        vyL  : (ivy, ivy-hshift_to_vxR)
        vyR  : (ivy, ivy+hshift_to_vxR)
        vxUL : (ivy, ivy-1-hshift_to_vxR)
        vxLL : (ivy, ivy+2-hshift_to_vxR)
        vxUR : (ivy, ivy-1)
        vxLR : (ivy, ivy+2)
        p'U  : (ivy, ivy+1)
        p'D  : (ivy, ivy+4)

    where ivy = ivx + 1

    For the Continuity equation coefficient matrix locations are defined as
    follows:

        vxL : (ipr, ipr-2-hshift_to_vxR)
        vxR : (ipr, ipr-2)
        vxU : (ipr, ipr-4)
        vyD : (ipr, ipr-1)



     j' 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27
        x  y  p  x  y  p  x  y  p  x  y  p  x  y  p  x  y  p  x  y  p  x  y  p  x  y  p
i' 1 x  c  LL L  D           UR    R  LR R
   2 y  UR c  U  LR D  D              R
   3 p  R  D  1
   4 x  U  UL    c  LL L  D           UR    R  LR R
   5 y     U     UR c  U  LR D  D              R
   6 p     U     R  D  0
   7 x           U  UL    c  LL L  D           UR    R  LR R
   8 y                       1
   9 p  L           U     R  D  0
  10 x  L                 U  UL    c  LL L  D           UR    R  LR R
  11 y  UL L     LL          U     UR c  U  LR D  D              R
  12 p           L           U     R  D  0
  13 x           L                 U  UL    c  LL L  D           UR    R  LR R
  14 y           UL L     LL          U     UR c  U  LR D  D              R
  15 p                    L           U     R  D  0
  16 x                    L                 U  UL    c  LL L  D           UR    R  LR R
  17 y                                                  1
  18 p                             L           U     R  D  0
  19 x                                                        1
  20 y                             UL L     LL          U     UR c  U  LR D  D
  21 p                                      L           U     R  D  0
  22 x                                                                 1
  23 y                                      UL L     LL          U     UR c  U  LR D  D
  24 p                                               L           U     R  D  0
  25 x                                                                          1
  26 y                                                                             1
  27 p                                                        L           U     R  D  0

        Figure 5: Sparse matrix for a basic grid with xnum = ynum = 4. "x",
        "y, and "p" denote rows and columns for unknowns for vx, vy and p',
        respectively. c denotes non-zero central coefficients. hshift_to_vxR
        is equal to 9. The number 1 in central node locations refers to where
        ghost boundary unknowns are defined for vx and vy unknowns. Zeros along
        the diagonal represent scaled pressure coefficients equal to zero.
        Blank spaces are zero-value coefficients for non-diagonal vx, vy and
        p' unknowns. Note that the p' central node has one non-zero
        coefficient at matrix location (2,2) where the pressure boundary
        condition is applied.

    Velocity Boundary Conditions
    ===========================================================================

    Velocity boundary condition are prescribed at the ghost velocity nodes
    using equations that relate velocity at ghost boundary nodes to velocity
    at unknowns. Coefficients for these equations have the following general
    form:

        vgn = C + D*v

    where vgn is the velocity component at the ghost node, v is the velocity
    unknown, C is prescribed velocity component and D is an integer. D can
    have a value of -1, 0 or 1 depending on the type of boundary condition.
    Boundary condition equation coefficients are stored in arrays that have
    the following nomenclature and structure:

        C = b<boundary_type><coordinate direction>[xnum+1 or ynum+1, 1]
        D = b<boundary_type><coordinate direction>[xnum+1 or ynum+1, 2]

    where <boundary_type> can be top, bottom, left or right, and
    <coordinate_direction> can be x or y. Note that the size of the array btopx
    is xnum+1 to account for xnum+1 staggered vy nodes in the x-direction even
    though only xnum staggered vx nodes are present. This allows the array to
    store coefficients for both vx and vy staggered grids. A similar
    explanation can be applied to the other arrays to account for array size.

    Boundary velocity equations using parameter storage arrays are as follows:

        Top boundary
        ------------
        Vx[1,j] = btopx[j,1] + btopx[j,2]*Vx[2,j]
        Vy[1,j] = btopy[j,1] + btopy[j,2]*Vy[2,j]

        Bottom boundary
        ---------------
        Vx[ynum+1,j] = bbottomx[j,1] + bbottomx[j,2]*Vx[ynum,j]
        Vy[ynum,j] = bbottomy[j,1] + bbottomy[j,2]*Vy[ynum-1,j]

        Left boundary
        --------------
        Vx[i,1] = bleftx[i,1] + bleftx[i,2]*Vx[i,2]
        Vy[i,1] = blefty[i,1] + blefty[i,2]*Vy[i,2]

        Right boundary
        ---------------
        Vx[i,xnum] = brightx[i,1] + brightx[i,2]*Vx[i,xnum-1]
        Vy[i,xnum+1] = brighty[i,1] + brighty[i,2]*Vy[i,xnum]

    No-slip boundary conditions are implemented for ghost nodes outside of the
    model domain by setting ghost node velocity equal to the opposite of the
    internal velocity making the velocity at the model boundary zero using
    D = -1. Velocity is set to zero using C=D=0 for ghost boundary nodes
    located along the boundary of the model domain. If no-slip conditions are
    applied along all boundaries, boundary equations would be as follows:

        Top boundary
        ------------
        Vx[1,j] = -Vx[2,j]
        Vy[1,j] = 0

        Bottom boundary
        --------------
        Vx[ynum+1,j] = -Vx[ynum,j]
        Vy[ynum,j] = 0

        Left boundary
        --------------
        Vx[i,1] = 0
        Vy[i,1] = -Vy[i,2]

        Right boundary
        ---------------
        Vx[i,xnum] = 0
        Vy[i,xnum+1] = -Vy[i,xnum]

    Free-slip boundary conditions with prescribed inflow/outflow are
    implemented by setting C=0 and D=1 for the boundary equation with
    components parallel to the boundary and D=0 and C equal to the
    inflow/outflow velocity for the boundary equation with components
    perpendicular to boundary. If free-slip conditions are applied along
    all boundaries, boundary equations would be as follows:

        Top boundary
        ------------
        Vx[1,j] = Vx[2,j]
        Vy[1,j] = C_in_out_y

        Bottom boundary
        ---------------
        Vx[ynum+1,j] = Vx[ynum,j]
        Vy[ynum,j] = C_in_out_y

        Left boundary
        --------------
        Vx[i,1] = C_in_out_x
        Vy[i,1] = Vy[i,2]

        Right boundary
        ---------------
        Vx[i,xnum] = C_in_out_x
        Vy[i,xnum+1] = Vy[i,xnum]

    Pressure boundary conditions
    ===========================================================================
    The following pressure boundary-condition options are available in this
    function:

        bpres_itype = 0: pressure is defined in basic grid cell located in the
                         upper left corner.

        bpres_itype = 1: pressure is defined along the top and bottom
                         boundaries with the top pressure set to the prescribed
                         pressure value and the bottom pressure set to zero.

        bpres_itype = 2: pressure is defined along the left and right
                         boundaries with the left pressure set to the
                         prescribed pressure value and the right pressure set
                         to zero.

    The value of the pressure boundary condition is controlled with the
    input parameter called pressure_bc.

    Internal prescribed velocity
    ===========================================================================

    Internal prescribed velocity is defined using the bintern_zone array and the
    bintern_velocity array:

        bintern_zone[1] : Horizontal position of vx nodes with prescribed velocity
                     [no such condition if negative]
        bintern_zone[2] : Min vertical position index i of zone with prescribed
                     velocity
        bintern_zone[3] : Max vertical position index i of zone with prescribed
                     velocity
        bintern_zone[4] : (not used)
        bintern_zone[5] : Horizontal position of vy nodes with prescribed velocity
                     [no such condition if negative]
        bintern_zone[6] : Min vertical position index j of zone with prescribed
                     velocity
        bintern_zone[7] : Max vertical position index j of zone with prescribed
                     velocity
        bintern_zone[8] : (not used)

        bintern_velocity[4] : Prescribed velocity in x-direction, m/s
        bintern_velocity[4] : Prescribed velocity in y-direction, m/s

"""
module StokesBuildManager

include("utils/BuildStructs.jl")
include("discretize_x_stokes_equation/StencilForVxStokesUnknown.jl")
include("discretize_y_stokes_equation/StencilForVyStokesUnknown.jl")
include("discretize_continuity_equation/StencilForPContinuityUnknown.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox: BuildSysTools
import EarthBox.BuildSysTools: SystemVectors
import EarthBox.Arrays: ArrayUtils
import .BuildStructs: CellIndices
import .BuildStructs: StokesBuildData
import .StencilForVxStokesUnknown
import .StencilForVyStokesUnknown
import .StencilForPContinuityUnknown

"""
    build_system_of_equations(model::ModelData)

This function loops over basic grid cells and calculates a discretized
system of equations using a conservative finite-difference approximation
on a staggered 2D grid with velocity ghost nodes (Figure 1). Only the
non-zero coefficients of matrix L are computed by this function. The matrix
indices of the non-zero elements, the basic grid indices Lii and Ljj
associated with each matrix element and the right-and-side vector R are
also computed.

The loop used in this function is over all solution nodes (See Figure 1
from stokes_build.py). Solution nodes are located in the upper
left-hand-corner of each basic grid cell.

# Arguments
- `model::ModelData`: Main model container object.

# Returns
- `system_vectors::SystemVectors`:
    - System vectors with non-zero matrix elements.
"""
function build_system_of_equations(model::ModelData)::SystemVectors
    build_data = StokesBuildData(model)
    xnum = build_data.grid.xnum
    ynum = build_data.grid.ynum
    inz = 1
    for i in 1:(ynum-1)
        for j in 1:(xnum-1)
            cell_indices = CellIndices(i, j, ynum)
            inz = apply_stencils(inz, cell_indices, build_data)
        end
    end
    # Subtract 1 from inz to account for the initial value of 1.
    return BuildSysTools.clean_non_zero_arrays!(inz-1, build_data.system_vectors)
end

""" Apply stencils to calculate non-zero matrix coefficients and rhs vector.

This function applies the stencils to the current cell and updates the
number of non-zero matrix elements.

# Arguments
- `inz::Int`: 
    - Number of non-zero matrix elements.

- `cell_indices::CellIndices`: 
    - Cell indices of the current cell.

- `build_data::StokesBuildData`: 
    - Build data.

# Returns
- `inz::Int`: 
    - Index of non-zero
      matrix arrays (Lii, Ljj, Li, Lj and Lv).

"""
function apply_stencils(
    inz::Int, 
    cell_indices::CellIndices, 
    build_data::StokesBuildData
)::Int
    inz = StencilForVxStokesUnknown.apply_stencil(inz, cell_indices, build_data)
    inz = StencilForVyStokesUnknown.apply_stencil(inz, cell_indices, build_data)
    inz = StencilForPContinuityUnknown.apply_stencil(inz, cell_indices, build_data)
    return inz
end

end # module
