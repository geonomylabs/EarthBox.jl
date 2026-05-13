""" Build compact system of equations for Stokes-continuity equations in 3D.

Description of the 3D Staggered Grid
================================================================================

The 3D staggered grid is composed of a basic grid and 4 staggered sub-grids.
The basic grid has dimensions (ynum, xnum, znum) and is indexed by (i, j, k)
where i is the vertical y-direction, j is the horizontal x-direction and k is
the horizontal z-direction. Basic-grid nodes sit at cell corners.

Because the grid cannot be drawn in a single ASCII figure, the geometry is
described with three orthogonal cross-sections. All three are drawn from
the same small illustrative grid: `ynum = 3`, `xnum = 4`, `znum = 3`, which
gives `(ynum-1) * (xnum-1) * (znum-1) = 12` interior cells and `N = 48`
unknowns. Every integer label that follows is derived from these sizes via
the cell-index formula in Eq. 15' below; in particular at slab `k = 1`:

    cell_index of (i, j, 1) = (j-1)*2 + i,    so cells 1..6 fill slab k=1
    cell_index of (i, j, 2) = 6 + (j-1)*2 + i, so cells 7..12 fill slab k=2

    ivx = 4*cell_index - 3,  ivy = ivx + 1,  ivz = ivx + 2,  ipr = ivx + 3

Throughout: "+" marks basic nodes; `vx`, `vy`, `vz` mark staggered velocity
nodes; `p'` marks scaled-pressure nodes. The `vx`/`vy`/`vz` token always
labels the velocity node itself — the basic-grid coordinate is given
separately as `i=…`, `j=…`, `k=…` row/column markers and as the digit in
the upper-left of each cell, which is the cell-index. Any digit appended
*to* a velocity token (e.g. `vy10`, `vx17`, `vz3`) is the matrix-row
index of that unknown.

The cross-section figures are 2D projections: the `vx`, `vy` and `p'`
nodes drawn in Cross-section A all actually sit at `z = k + 0.5` (mid-
slab), half a basic-slab away from the `+` corners which lie on slabs
`z = k` and `z = k+1`. Only `vz` (visible in Cross-section B) is co-
planar with a basic slab. Analogous projection offsets apply to Cross-
sections B and C. Ghost-node columns/rows outside the basic grid are
drawn at the same horizontal positions as the labeled internal nodes so
that each ghost is visually unambiguous as a node.

::

    Cross-section A — xy-plane through slab k = 1

                      j=1             j=2             j=3             j=4
                            xstp(1)         xstp(2)         xstp(3)
                       |...............|...............|...............|

                    xstpc(1)        xstpc(2)        xstpc(3)        xstpc(4)
               |...............|...............|...............|...............|

                       vx              vx              vx              vx                 ......  <- top y-ghost vx
       
        
        i=1    vy       +------vy-------+------vy-------+------vy-------+      vy ......  ystpc(1)  (top y-boundary)
                       |1              |3              |5              |
                       |               |               |               |
                       vx     p'4     vx1     p'12   vx9     p'20    vx17         ystp(1) ......
                       |               |               |               |
                       |               |               |               |
        i=2    vy       +-----vy2-------+-----vy10------+-----vy18------+      vy ......  ystpc(2)
                       |2              |4              |6              |
                       |               |               |               |
                       vx     p'8     vx5     p'16   vx13    p'24    vx21         ystp(2) ......
                       |               |               |               |
                       |               |               |               |
        i=3    vy       +-----vy6-------+-----vy14------+-----vy22------+      vy ......  ystpc(3)  (bottom y-boundary)


                       vx              vx              vx              vx                 ......  <- bottom y-ghost vx

    Row "i=1" is the top y-boundary: the `vy` nodes *on* this row are
    ghost unknowns (no matrix index) because Vy[1, j, k] maps to
    Vy[2, j, k] via `btopy`. Row "i=3" is the bottom y-boundary: the
    `vy` nodes there *do* carry matrix indices (`vy6`, `vy14`, `vy22`)
    because Vy[ynum, j, k] is the vy of cells at `i = ynum-1` and those
    cells' rows are filled by the y-Stokes BoundaryCells path.

    The `vy` tokens *outside* the leftmost and rightmost `+` on every
    basic-row line are the x-direction Vy ghosts (`Vy[i, 1, k]` on the
    left, `Vy[i, xnum+1, k]` on the right). Vy is staggered in `j`, so
    these positions sit outside the basic grid in x and are unlabeled
    ghost unknowns mapped to the adjacent interior `vy` via `bleftx` /
    `brightx`. The leftmost `vx` column (at `j=1`) is the left
    x-boundary vx (no matrix index); `vx17` on the right (at `j=xnum`)
    is the right x-boundary vx, carried by cells at `j = xnum-1`. The
    top and bottom `vx` ghost rows are drawn at the same `j` positions
    as the interior `vx` columns.

::

    Cross-section B — xz-plane through cell strip i = 1
    (slice between basic rows 1 and 2 of grid in Cross-section A)

                            j=1             j=2             j=3             j=4
                                  xstp(1)         xstp(2)         xstp(3)
                             |...............|...............|...............|
                          xstpc(1)        xstpc(2)        xstpc(3)        xstpc(4)
                     |...............|...............|...............|...............|

                             vx              vx              vx              vx                ......  <- front z-ghost vx


        k=1         vz       +------vz-------+------vz-------+------vz-------+      vz ......  zstpc(1)  (front boundary)
                             |1              |3              |5              |
                             |               |               |               |
                             vx     p'4     vx1     p'12   vx9     p'20    vx17        zstp(1) ......
                             |               |               |               |
                             |               |               |               |
        k=2         vz       +-----vz3-------+-----vz11------+-----vz19------+      vz ......  zstpc(2)
                             |7              |9              |11             |
                             |               |               |               |
                             vx     p'28   vx25     p'36   vx33    p'44    vx41        zstp(2) ......
                             |               |               |               |
                             |               |               |               |
        k=3         vz       +-----vz27------+-----vz35------+-----vz43------+      vz ......  zstpc(3)  (back boundary)


                             vx              vx              vx              vx                ......  <- back z-ghost vx

    Slab `k=1` is the front z-boundary: vz nodes *on* this slab are
    ghosts (Vz[i, j, 1] = bfrontz + bfrontz·Vz[i, j, 2]). Slab `k=3`
    carries matrix indices: under the back-face convention `vz of cell
    (i, j, k) = Vz[i+1, j+1, k+1]`, so slab 3 holds vz of cells with
    `k = znum-1 = 2`, e.g. `vz27 = ivz of cell (1, 1, 2)`. Slab `k=2`
    holds vz of cells with `k = 1`, so `vz3 = ivz of cell (1, 1, 1)`.

    The `vz` tokens *outside* the leftmost and rightmost `+` on every
    basic-slab line are the x-direction Vz ghosts (`Vz[i, 1, k]` on the
    left, `Vz[i, xnum+1, k]` on the right). Vz is staggered in `j`, so
    these positions sit outside the basic grid in x and are unlabeled
    ghost unknowns mapped to the adjacent interior `vz` via `bleftx` /
    `brightx`. The interior strip between slabs 1 and 2 holds vx and p'
    for cells `(1, *, 1)`; the strip between slabs 2 and 3 holds them
    for cells `(1, *, 2)`.

::

    Cross-section C — yz-plane through cell strip j = 1
    (slice between basic cols 1 and 2 of the same grid)

                          (front)                          (back)

                            k=1             k=2             k=3
                                   zstp(1)         zstp(2)
                             |...............|...............|

                             vz              vz              vz                     ......  <- top y-ghost vz

        i=1         vy       +------vy-------+------vy-------+      vy ..... ystpc(1)
                             |1              |7              |
                             |               |               |
                             vz     p'4     vz3     p'28   vz27        ystp(1) ......
                             |               |               |
                             |               |               |
        i=2         vy       +-----vy2-------+-----vy26------+      vy ...... ystpc(2)
                             |2              |8              |
                             |               |               |
                             vz     p'8     vz7     p'32   vz31        ystp(2) ......
                             |               |               |
                             |               |               |
        i=3         vy       +-----vy6-------+-----vy30------+      vy ...... ystpc(3)

                             vz              vz              vz                     ......  <- bottom y-ghost vz


    Same boundary rules as before: row `i=1` = top y-boundary (vy on
    this row are ghosts); row `i=3` = bottom y-boundary (vy carries
    matrix indices of cells at `i = ynum-1`). Column `k=1` (slab 1) is
    the front z-boundary, so the vz nodes there are ghosts; column
    `k=3` (slab 3) is the back z-boundary, so vz there carries matrix
    indices of cells at `k = znum-1`.

    The `vy` tokens *outside* the leftmost and rightmost `+` on every
    basic-row line are the z-direction Vy ghosts (`Vy[i, j, 1]` on the
    left, `Vy[i, j, znum+1]` on the right). Vy is staggered in `k`, so
    these positions sit outside the basic grid in z and are unlabeled
    ghost unknowns mapped to the adjacent interior `vy` via `bfronty` /
    `bbacky`. The top and bottom `vz` ghost rows are similarly the
    y-direction Vz ghosts (`Vz[1, j, k]` and `Vz[ynum+1, j, k]`),
    aligned with the interior `vz` columns.

    The interior strip between `k=1` and `k=2` holds vy and p' for
    cells `(*, 1, 1)`; the strip between `k=2` and `k=3` holds them
    for cells `(*, 1, 2)`.

Figure 1 (A, B, C): three orthogonal cross-sections of the 3D staggered
grid. Sub-grid array indices match the matrix-row indices listed in Eqs.
15'–18 below: e.g. in Cross-section A the label `vx9` is the vx unknown
of cell (1, 2, 1), located at matrix row `ivx = 4·3 − 3 = 9`. The same
cell shows `vy10` (`ivy = 10`) on its bottom face, `p'12` (`ipr = 12`)
at its centre, and `vz11` (`ivz = 11`) on its back face (visible in
Cross-section B).

The four staggered sub-grids are:

    (1) Vx of dimensions (ynum+1, xnum,   znum+1) — vx is located at the
        center of each x-perpendicular face of a basic grid cell. The vx of
        cell (i, j, k) is Vx[i+1, j+1, k+1].
    (2) Vy of dimensions (ynum,   xnum+1, znum+1) — vy is located at the
        center of each y-perpendicular face. The vy of cell (i, j, k) is
        Vy[i+1, j+1, k+1] (the bottom face of the cell).
    (3) Vz of dimensions (ynum+1, xnum+1, znum)   — vz is located at the
        center of each z-perpendicular face. The vz of cell (i, j, k) is
        Vz[i+1, j+1, k+1] (the back face of the cell).
    (4) p' of dimensions (ynum-1, xnum-1, znum-1) — p' lives at the centers
        of the interior basic grid cells.

Velocity nodes that lie on or outside the model boundary are ghost boundary
nodes used to define boundary conditions. Numbered ghost boundary nodes are
ghost unknowns used as placeholders in the sparse matrix.

Grid spacings along each axis are stored in arrays:

    xstp(xnum-1), ystp(ynum-1), zstp(znum-1)

For staggered velocity faces the spacings centered on basic nodes are:

    xstpc(xnum),  ystpc(ynum),  zstpc(znum)

The following arrays store fields defined on the basic grid:

    rho(ynum, xnum, znum)      : density (kg/m^3)
    tk(ynum, xnum, znum)       : temperature (K)
    kt(ynum, xnum, znum)       : thermal conductivity (W/m/K)
    hr(ynum, xnum, znum)       : radiogenic heat production (W/m^3)

Three shear viscosities are required, one per shear-stress component, each
defined on the basic edge family where that stress lives:

    etas_xy(ynum, xnum,   znum-1) : viscosity for σxy on z-parallel edges
    etas_xz(ynum-1, xnum, znum)   : viscosity for σxz on y-parallel edges
    etas_yz(ynum,   xnum-1, znum) : viscosity for σyz on x-parallel edges

The 2D `etas` field corresponds to `etas_xy` in this 3D notation.

Normal viscosity and stress are defined at pressure nodes:

    etan(ynum-1, xnum-1, znum-1) : viscosity for normal stress (Pa.s)
    sxx, syy, szz                : normal stresses (Pa) at pressure nodes
    p(ynum-1, xnum-1, znum-1)    : pressure (Pa)
    p'(ynum-1, xnum-1, znum-1)   : scaled pressure

Velocity arrays:

    Vx(ynum+1, xnum,   znum+1) : x-component of velocity (m/s)
    Vy(ynum,   xnum+1, znum+1) : y-component of velocity (m/s)
    Vz(ynum+1, xnum+1, znum)   : z-component of velocity (m/s)

Unscaled pressure p is related to scaled pressure p' via:

    p = pscale*p'                                                    (Eq. 1)

with the 3D scaling

    pscale = 3.0*etan[1,1,1] / (xstpavr + ystpavr + zstpavr)        (Eq. 2')

where etan[1,1,1] is the normal viscosity at the first interior pressure
node and xstpavr, ystpavr, zstpavr are the average basic-grid spacings.

Discretizing Governing Equations Using Finite Difference Stencils
================================================================================

Coefficients and right-hand-side terms are generated by looping over basic
grid cells and discretizing the governing equations associated with the
unknowns vx, vy, vz and p' using the relevant stencil. Each basic cell now
contributes four equations (one per unknown), in contrast to the 2D case
where each cell contributed three.

x-Stokes stencil applied for each vx unknown
--------------------------------------------------------------------------------

The x-Stokes stencil is depicted in two coupled cross-sections that share
the central node vxC. Cross-section A is the xy-plane through vxC and is
structurally identical to Figure 2 in the 2D docstring. Cross-section B is
the xz-plane through vxC and carries the new ∂σxz/∂z contribution.

::

    x-Stokes stencil — Section A (xy-plane through vxC, matches 2D Fig. 2)
    +------------------------+------------------------+
    |                        |                        |
    |                    Vx(i,j+1,k+1)                |
    |                        vxU                      |
    |                        |                        |
    BN(i,j,k)---Vy(i,j+1,k+1)---etas_xy_U---Vy(i,j+2,k+1)----ystpc(i)
    |          vyUL          etas_xy_U          vyUR          dyU
    |                        rhoU                     |
    |       P'(i,j,k)        |        P'(i,j+1,k)    |
    Vx(i+1,j,k+1)   p'L  Vx(i+1,j+1,k+1)   p'R   Vx(i+1,j+2,k+1)--ystp(i)
    |  vxL       etanL       vxC       etanR       vxR         dyC
    |                        |                        |
    +-----Vy(i+1,j+1,k+1)--etas_xy_D--Vy(i+1,j+2,k+1)----------ystpc(i+1)
    |          vyLL          etas_xy_D         vyLR          dyD
    |                    Vx(i+2,j+1,k+1)              |
    |                        vxD                      |
    +------------------------+------------------------+

::

    x-Stokes stencil — Section B (xz-plane through vxC, NEW for 3D)
    +------------------------+------------------------+
    |                        |                        |
    |                    Vx(i+1,j+1,k)                |
    |                        vxF                      |
    |                        |                        |
    BN(i,j,k)----Vz(i+1,j+1,k)----etas_xz_F----Vz(i+1,j+2,k)----zstpc(k)
    |         vzFL          etas_xz_F          vzFR          dzF
    |                        |                        |
    Vx(i+1,j,k+1)            vxC            Vx(i+1,j+2,k+1)----- xstpc(j+1)
    |                        |                         |
    |                        |                        |
    +-----Vz(i+1,j+1,k+1)--etas_xz_B--Vz(i+1,j+2,k+1)----------zstpc(k+1)
    |         vzBL          etas_xz_B          vzBR          dzB
    |                        |                        |
    |                    Vx(i+1,j+1,k+2)               |
    |                        vxB                       |
    +------------------------+------------------------+

The x-Stokes equation expresses force balance per unit volume at vxC:

    dσxx/dx + dσxy/dy + dσxz/dz - dP/dx = RX                       (Eq. 3')

with RX the right-hand-side. Using the two cross-sections above the
discretized form is:

      [ 2*etanR*(vxR - vxC)/dxR - 2*etanL*(vxC - vxL)/dxL ] / dxC
    + { etas_xy_D*[ (vxD - vxC)/dyD + (vyLR - vyLL)/dxC ]
       - etas_xy_U*[ (vxC - vxU)/dyU + (vyUR - vyUL)/dxC ] } / dyC
    + { etas_xz_B*[ (vxB - vxC)/dzB + (vzBR - vzBL)/dxC ]
       - etas_xz_F*[ (vxC - vxF)/dzF + (vzFR - vzFL)/dxC ] } / dzC
    - pscale*(p'R - p'L)/dxC = -gravity_x*(rhoU+rhoD)/2            (Eq. 4')

Rearranging gives the coefficient list (suffix legend below):

    - (   2*etanL/dxL/dxC + 2*etanR/dxR/dxC
        + etas_xy_U/dyU/dyC + etas_xy_D/dyD/dyC
        + etas_xz_F/dzF/dzC + etas_xz_B/dzB/dzC ) * vxC
    + 2*etanL/dxL/dxC * vxL
    + 2*etanR/dxR/dxC * vxR
    + etas_xy_U/dyU/dyC * vxU
    + etas_xy_D/dyD/dyC * vxD
    + etas_xz_F/dzF/dzC * vxF
    + etas_xz_B/dzB/dzC * vxB
    + etas_xy_U/dxC/dyC * vyUL  - etas_xy_U/dxC/dyC * vyUR
    - etas_xy_D/dxC/dyC * vyLL  + etas_xy_D/dxC/dyC * vyLR
    + etas_xz_F/dxC/dzC * vzFL  - etas_xz_F/dxC/dzC * vzFR
    - etas_xz_B/dxC/dzC * vzBL  + etas_xz_B/dxC/dzC * vzBR
    - pscale/dxC * p'R
    + pscale/dxC * p'L
    = -gravity_x*(rhoU+rhoD)/2                                      (Eq. 7')

Suffix legend (additions to 2D):

    C = central, L/R = -/+ x, U/D = -/+ y, F/B = -/+ z.
    Compound suffixes follow the same convention, e.g. UFL = upper-front-left.

The x-Stokes stencil is not applied along basic cells lying on the right
boundary (j = xnum-1), where vxC is prescribed. In that case the matrix
element (ivx, ivx) is set to a Kbond-scaled 1 and the right-hand-side is
set to the prescribed velocity, as in 2D. Ghost-boundary nodes appearing
in the stencil cause some unknowns to become constants, triggering the
same regrouping and rhs adjustment described in the 2D docstring.

y-Stokes stencil applied for each vy unknown
--------------------------------------------------------------------------------

The y-Stokes stencil is depicted in two coupled cross-sections that share
the central node vyC. Cross-section A is the xy-plane through vyC and is
structurally identical to Figure 3 in the 2D docstring. Cross-section B
is the yz-plane through vyC and carries the new ∂σyz/∂z contribution with
shear viscosity etas_yz on x-parallel edges.

::

    y-Stokes stencil — Section A (xy-plane through vyC, matches 2D Fig. 3)
    +-------------------BN(i,j,k)---Vy(i,j+1,k+1)----------+----------------------+----
    |                         |          vyU               |                      |
    |                         |                            |                      |
    |                         |       P'(i,j,k)            |                      |
    |                    Vx(i+1,j,k+1)    p'U        Vx(i+1,j+1,k+1)              |ystp(i)----
    |                         vxUL    etan(i,j,k)        vxUR                     | dyU
    |                         |          etanU             |                      |
    |                         |                            |                      |
    +----Vy(i+1,j,k+1)---etas_xy_L---Vy(i+1,j+1,k+1)---etas_xy_R---Vy(i+1,j+2,k+1)+--ystpc(i+1)
    |       vyL         etas_xy_L         vyC          etas_xy_R       vyR        | dyC
    |                         rhoL                         rhoR                   |
    |                         |       P'(i+1,j,k)          |                      |
    |                         |          p'D               |                      |
    |                    Vx(i+2,j,k+1)  etan(i+1,j,k)  Vx(i+2,j+1,k+1)            |ystp(i+1)---
    |                         vxLL        etanD          vxLR                     | dyD
    |                         |                            |                      |
    +-------------------------+-------Vy(i+2,j+1,k+1)-----+----------------------+----
                                         vyD

           |          xstpc(j)        |          xstpc(j+1)         |
                         dxL                          dxR
                          |          xstp(j)          |
                                       dxC

::

    y-Stokes stencil — Section B (yz-plane through vyC, NEW for 3D)
    +-----------------BN(i,j+1,k)---Vy(i,j+1,k+1)----------+----------------------+----
    |                         |          vyU               |                      |
    |                         |                            |                      |
    |                         |       P'(i,j,k)            |                      |
    |                    Vz(i+1,j+1,k)    p'U         Vz(i+1,j+1,k+1)             |ystp(i)----
    |                         vzFU    etan(i,j,k)        vzBU                     | dyU
    |                         |          etanU             |                      |
    |                         |                            |                      |
    +----Vy(i+1,j+1,k)---etas_yz_F---Vy(i+1,j+1,k+1)---etas_yz_B---Vy(i+1,j+1,k+2)+--ystpc(i+1)
    |       vyF         etas_yz_F         vyC          etas_yz_B       vyB        | dyC
    |                         |       P'(i+1,j,k)          |                      |
    |                         |          p'D               |                      |
    |                    Vz(i+2,j+1,k)  etan(i+1,j,k)  Vz(i+2,j+1,k+1)            |ystp(i+1)---
    |                         vzFD        etanD          vzBD                     | dyD
    |                         |                            |                      |
    +-------------------------+-------Vy(i+2,j+1,k+1)-----+----------------------+----
                                         vyD

           |          zstpc(k)        |          zstpc(k+1)         |
                         dzF                          dzB
                          |          zstp(k)          |
                                       dzC

The y-Stokes equation expresses force balance per unit volume at vyC:

    dσyy/dy + dσyx/dx + dσyz/dz - dP/dy = RY                       (Eq. 8')

with RY the right-hand-side. Using the two cross-sections above the
discretized form is:

      [ 2*etanD*(vyD - vyC)/dyD - 2*etanU*(vyC - vyU)/dyU ] / dyC
    + { etas_xy_R*[ (vyR - vyC)/dxR + (vxLR - vxUR)/dyC ]
       - etas_xy_L*[ (vyC - vyL)/dxL + (vxLL - vxUL)/dyC ] } / dxC
    + { etas_yz_B*[ (vyB - vyC)/dzB + (vzBD - vzBU)/dyC ]
       - etas_yz_F*[ (vyC - vyF)/dzF + (vzFD - vzFU)/dyC ] } / dzC
    - pscale*(p'D - p'U)/dyC = -gravity_y * <rho>                  (Eq. 9')

where `<rho>` is the average of the four basic-grid densities at the
corners of the vyC face — rho[i+1,j,k], rho[i+1,j+1,k], rho[i+1,j,k+1]
and rho[i+1,j+1,k+1] — i.e. an x- and z-averaging of `rhoL`/`rhoR` from
Section A across the two slabs in Section B.

Rearranging gives the coefficient list (suffix legend identical to the
x-Stokes legend above):

    - (   2*etanU/dyU/dyC + 2*etanD/dyD/dyC
        + etas_xy_L/dxL/dxC + etas_xy_R/dxR/dxC
        + etas_yz_F/dzF/dzC + etas_yz_B/dzB/dzC ) * vyC
    + 2*etanU/dyU/dyC * vyU
    + 2*etanD/dyD/dyC * vyD
    + etas_xy_L/dxL/dxC * vyL
    + etas_xy_R/dxR/dxC * vyR
    + etas_yz_F/dzF/dzC * vyF
    + etas_yz_B/dzB/dzC * vyB
    + etas_xy_L/dxC/dyC * vxUL  - etas_xy_R/dxC/dyC * vxUR
    - etas_xy_L/dxC/dyC * vxLL  + etas_xy_R/dxC/dyC * vxLR
    + etas_yz_F/dyC/dzC * vzFU  - etas_yz_B/dyC/dzC * vzBU
    - etas_yz_F/dyC/dzC * vzFD  + etas_yz_B/dyC/dzC * vzBD
    - pscale/dyC * p'D
    + pscale/dyC * p'U
    = -gravity_y * <rho>                                            (Eq. 10')

Gravity acts along y, so the y-Stokes stencil is the only one that
carries an interface-stabilization correction (Duretz et al., 2011). The
correction adds `-drhody_gy_dt` to the diagonal vy coefficient,
`-drhodx_gy_dt/4` to each of the four vx edge neighbours (vxUL, vxUR,
vxLL, vxLR), and `-drhodz_gy_dt/4` to each of the four vz edge
neighbours (vzFU, vzFD, vzBU, vzBD).

The y-Stokes stencil is not applied along basic cells lying on the
bottom boundary (i = ynum-1), where vyC is prescribed. In that case the
matrix element (ivy, ivy) is set to a Kbond-scaled 1 and the right-hand-
side is set to the prescribed velocity, as in 2D. Ghost-boundary nodes
appearing in the stencil cause some unknowns to become constants,
triggering the same regrouping and rhs adjustment described in the 2D
docstring.

z-Stokes stencil applied for each vz unknown  (NEW)
--------------------------------------------------------------------------------

Symmetric to x- and y-Stokes. Section A is the xz-plane through vzC and
carries the ∂σxz/∂x contribution (via etas_xz). Section B is the yz-plane
through vzC and carries the ∂σyz/∂y contribution (via etas_yz). The
z-Stokes equation is:

    dσzz/dz + dσzx/dx + dσzy/dy - dP/dz = RZ                       (Eq. 8'')

with discretized form analogous to Eq. 4' but centered on vzC. The
coefficient list has the same form as Eq. 7' / Eq. 10 with axis labels
rotated so the dominant normal-stress neighbours are vzF, vzB and the
pressure pair is p'F, p'B.

The z-Stokes stencil is not applied along basic cells lying on the back
boundary (k = znum-1), where vzC is prescribed; the same Kbond-scaled
prescribed-velocity treatment is used.

Continuity stencil applied for each p' unknown
--------------------------------------------------------------------------------

::

    Continuity stencil — Section A (xy-plane through p'C, matches 2D Fig. 4)

                            BN(i,j,k)------Vy(i,j+1,k+1)------+
                                |             vyU             |
                                |                             |
                                |                             |
                                |                             |
                            Vx(i+1,j,k+1)   P'(i,j,k)   Vx(i+1,j+1,k+1)   ystp(i)
                                vxL            p'C            vxR           dy
                                |                             |
                                |                             |
                                |                             |
                                +-----Vy(i+1,j+1,k+1)---------+
                                            vyD

                                |........xstp(j)............|
                                            dx

::

    Continuity stencil — Section B (xz-plane through p'C, NEW for 3D)

                            BN(i,j,k)------Vz(i+1,j+1,k)------+
                                |             vzF             |
                                |                             |
                                |                             |
                                |                             |
                            Vx(i+1,j,k+1)   P'(i,j,k)   Vx(i+1,j+1,k+1)   zstp(k)
                                vxL            p'C            vxR           dz
                                |                             |
                                |                             |
                                |                             |
                                +-----Vz(i+1,j+1,k+1)---------+
                                            vzB

                                |........xstp(j)............|
                                            dx

    Figure 4': continuity equation stencil in 3D. Section A is structurally
    identical to the 2D continuity figure and provides the (vxR − vxL)/dx
    and (vyD − vyU)/dy terms. Section B is the new xz cut through p'C and
    provides the (vzB − vzF)/dz term. (A symmetric yz cut through p'C would
    provide the same dvy/dy term as Section A and the same dvz/dz term as
    Section B, so it is omitted.)

The 3D continuity equation is:

    dvx/dx + dvy/dy + dvz/dz = RC                                  (Eq. 11')

Discretized:

    pscale*[ (vxR - vxL)/dx + (vyD - vyU)/dy + (vzB - vzF)/dz ] = RC*pscale

Rearranged coefficient list:

      pscale/dx * vxR  - pscale/dx * vxL
    + pscale/dy * vyD  - pscale/dy * vyU
    + pscale/dz * vzB  - pscale/dz * vzF
    = RC*pscale

Building the System of Equations
================================================================================

A sparse system

    LX = R                                                          (Eq. 13)

is built with size N x N, where

    N = (xnum-1) * (ynum-1) * (znum-1) * 4                         (Eq. 14')

reflecting four unknowns per interior basic cell (vx, vy, vz, p'). The
central node for a cell (i, j, k) appears in the large matrix at indices:

    cell_index = ((k-1)*(xnum-1) + (j-1))*(ynum-1) + i             (Eq. 15')

    ivx = cell_index*4 - 3                                          (Eq. 16'a)
    ivy = ivx + 1                                                   (Eq. 16'b)
    ivz = ivx + 2                                                   (Eq. 16'c)
    ipr = ivx + 3                                                   (Eq. 16'd)

Two strides are used to address off-diagonal neighbours instead of one:

    hshift_to_vxR = (ynum-1) * 4                                   (Eq. 17')
    dshift_to_vxF = (xnum-1) * (ynum-1) * 4                        (Eq. 18 NEW)

`hshift_to_vxR` is the row stride from a cell to its j+1 neighbour (the
2D shift, with 4 replacing 3). `dshift_to_vxF` is the new row stride from
a cell to its k+1 neighbour.

x-Stokes coefficient matrix locations:

    vxC  : (ivx, ivx)
    vxL  : (ivx, ivx - hshift_to_vxR)
    vxR  : (ivx, ivx + hshift_to_vxR)
    vxU  : (ivx, ivx - 4)
    vxD  : (ivx, ivx + 4)
    vxF  : (ivx, ivx - dshift_to_vxF)
    vxB  : (ivx, ivx + dshift_to_vxF)
    vyUL : (ivx, ivx - 3)
    vyUR : (ivx, ivx - 3 + hshift_to_vxR)
    vyLL : (ivx, ivx + 1)
    vyLR : (ivx, ivx + 1 + hshift_to_vxR)
    vzFL : (ivx, ivx + 2 - dshift_to_vxF)
    vzFR : (ivx, ivx + 2 + hshift_to_vxR - dshift_to_vxF)
    vzBL : (ivx, ivx + 2)
    vzBR : (ivx, ivx + 2 + hshift_to_vxR)
    p'L  : (ivx, ivx + 3)
    p'R  : (ivx, ivx + 3 + hshift_to_vxR)

y-Stokes and z-Stokes coefficient locations follow the same pattern with
axis labels rotated; ivy = ivx + 1 and ivz = ivx + 2 in all expressions.

Continuity coefficient locations:

    vxL : (ipr, ipr - 3 - hshift_to_vxR)
    vxR : (ipr, ipr - 3)
    vyU : (ipr, ipr - 6)
    vyD : (ipr, ipr - 2)
    vzF : (ipr, ipr - 1 - dshift_to_vxF)
    vzB : (ipr, ipr - 1)

The sparse pattern is too wide to depict in ASCII for any meaningful 3D
case; conceptually it has the same band structure as Figure 5 of the 2D
docstring, with additional bands offset by ±dshift_to_vxF carrying the
k±1 neighbour contributions. Central coefficients lie at (i', i') in
blocks of 4 (vx, vy, vz, p').

Velocity Boundary Conditions
================================================================================

Velocity boundary conditions are prescribed at ghost velocity nodes via
equations of the form

    vgn = C + D*v

where vgn is the velocity component at the ghost node, v is the velocity
unknown, C is a prescribed component and D ∈ {-1, 0, 1}. Coefficients are
stored in arrays named

    C = b<boundary_type><coordinate_direction>[ ..., 1 ]
    D = b<boundary_type><coordinate_direction>[ ..., 2 ]

where <boundary_type> ∈ {top, bottom, left, right, front, back} and
<coordinate_direction> ∈ {x, y, z}. Each boundary array is a 3D structure:
two indices address a node on the 2D face plus a final size-2 axis for
the (C, D) pair.

Per-face boundary equations:

    Top boundary    (i = 1 ghost row)
        Vx[1,j,k] = btopx[j,k,1] + btopx[j,k,2] * Vx[2,j,k]
        Vy[1,j,k] = btopy[j,k,1] + btopy[j,k,2] * Vy[2,j,k]
        Vz[1,j,k] = btopz[j,k,1] + btopz[j,k,2] * Vz[2,j,k]

    Bottom boundary (i = ynum+1 ghost row for Vx and Vz, i = ynum for Vy)
        Vx[ynum+1,j,k] = bbottomx[j,k,1] + bbottomx[j,k,2] * Vx[ynum,j,k]
        Vy[ynum,j,k]   = bbottomy[j,k,1] + bbottomy[j,k,2] * Vy[ynum-1,j,k]
        Vz[ynum+1,j,k] = bbottomz[j,k,1] + bbottomz[j,k,2] * Vz[ynum,j,k]

    Left boundary   (j = 1 ghost column)
        Vx[i,1,k] = bleftx[i,k,1] + bleftx[i,k,2] * Vx[i,2,k]
        Vy[i,1,k] = blefty[i,k,1] + blefty[i,k,2] * Vy[i,2,k]
        Vz[i,1,k] = bleftz[i,k,1] + bleftz[i,k,2] * Vz[i,2,k]

    Right boundary  (j = xnum for Vx, j = xnum+1 ghost column for Vy/Vz)
        Vx[i,xnum,k]   = brightx[i,k,1] + brightx[i,k,2] * Vx[i,xnum-1,k]
        Vy[i,xnum+1,k] = brighty[i,k,1] + brighty[i,k,2] * Vy[i,xnum,k]
        Vz[i,xnum+1,k] = brightz[i,k,1] + brightz[i,k,2] * Vz[i,xnum,k]

    Front boundary  (k = 1 ghost slab)   [NEW]
        Vx[i,j,1] = bfrontx[i,j,1] + bfrontx[i,j,2] * Vx[i,j,2]
        Vy[i,j,1] = bfronty[i,j,1] + bfronty[i,j,2] * Vy[i,j,2]
        Vz[i,j,1] = bfrontz[i,j,1] + bfrontz[i,j,2] * Vz[i,j,2]

    Back boundary   (k = znum for Vz, k = znum+1 ghost slab for Vx/Vy) [NEW]
        Vx[i,j,znum+1] = bbackx[i,j,1] + bbackx[i,j,2] * Vx[i,j,znum]
        Vy[i,j,znum+1] = bbacky[i,j,1] + bbacky[i,j,2] * Vy[i,j,znum]
        Vz[i,j,znum]   = bbackz[i,j,1] + bbackz[i,j,2] * Vz[i,j,znum-1]

No-slip and free-slip cases reduce as in 2D:

  - No-slip:    components parallel to the face use D = -1 (mirror, zero
                velocity on the boundary); the component perpendicular to
                the face uses C = D = 0 (zero on the boundary itself).
  - Free-slip: components parallel to the face use D = +1 (zero gradient);
                the perpendicular component uses D = 0 and C equal to any
                prescribed inflow/outflow.

Pressure Boundary Conditions
================================================================================

Four pressure boundary-condition modes are available:

    bpres_itype = 0: pressure is defined in the basic grid cell located at
                     the upper-left-front corner (i = j = k = 1).

    bpres_itype = 1: pressure is defined along the top and bottom (y)
                     boundaries; top pressure is set to the prescribed
                     value and bottom pressure is set to zero.

    bpres_itype = 2: pressure is defined along the left and right (x)
                     boundaries; left pressure is set to the prescribed
                     value and right pressure is set to zero.

    bpres_itype = 3 (NEW): pressure is defined along the front and back
                     (z) boundaries; front pressure is set to the
                     prescribed value and back pressure is set to zero.

The value applied is controlled by the input parameter `pressure_bc`.

Internal Prescribed Velocity
================================================================================

Internal prescribed velocity is defined with the `bintern_zone` array and
the `bintern_velocity` array using a 3D layout that extends the 2D
convention. Each prescribed-velocity zone uses six consecutive entries:

  bintern_zone[1..6]  : vx zone — [j_col, min_i, max_i, min_k, max_k, unused]
                        (no prescription if bintern_zone[1] is negative)
  bintern_zone[7..12] : vy zone — [i_row, min_j, max_j, min_k, max_k, unused]
                        (no prescription if bintern_zone[7] is negative)
  bintern_zone[13..18]: vz zone — [k_slab, min_j, max_j, min_i, max_i, unused]
                        (no prescription if bintern_zone[13] is negative)

  bintern_velocity[1] : prescribed vx velocity (m/s)
  bintern_velocity[2] : prescribed vy velocity (m/s)
  bintern_velocity[3] : prescribed vz velocity (m/s)
"""
module StokesBuildManager3D

include("utils/BuildStructs.jl")
include("utils/Predicates3D.jl")
include("utils/DebugStencils.jl")
include("discretize_x_stokes_equation/StencilForVxStokesUnknown.jl")
include("discretize_y_stokes_equation/StencilForVyStokesUnknown.jl")
include("discretize_z_stokes_equation/StencilForVzStokesUnknown.jl")
include("discretize_continuity_equation/StencilForPContinuityUnknown.jl")

import EarthBox.ModelDataContainer: ModelData
import EarthBox: BuildSysTools
import EarthBox.BuildSysTools: SystemVectors
import EarthBox.Arrays: ArrayUtils
import .BuildStructs: CellIndices
import .BuildStructs: StokesBuildData
import .StencilForVxStokesUnknown
import .StencilForVyStokesUnknown
import .StencilForVzStokesUnknown
import .StencilForPContinuityUnknown

"""
    build_system_of_equations(model::ModelData)

This function loops over basic grid cells and calculates a discretized
system of equations using a conservative finite-difference approximation
on a staggered 3D grid with velocity ghost nodes (Figure 1). Only the
non-zero coefficients of matrix L are computed by this function. The
matrix indices of the non-zero elements, the basic grid indices Lii and
Ljj associated with each matrix element and the right-hand-side vector R
are also computed.

The loop visits every basic grid cell (i, j, k) for i in 1..ynum-1, j in
1..xnum-1, k in 1..znum-1. The loop order is `i` innermost, then `j`,
then `k`, matching the column-major basic-grid layout and the cell-index
formula in Eq. 15'.

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
    znum = build_data.grid.znum
    inz = 1
    @inbounds for k in 1:(znum-1)
        for j in 1:(xnum-1)
            for i in 1:(ynum-1)
                cell_indices = CellIndices(i, j, k, ynum, xnum)
                inz = apply_stencils(inz, cell_indices, build_data)
            end
        end
    end
    # Subtract 1 from inz to account for the initial value of 1.
    return BuildSysTools.clean_non_zero_arrays!(inz-1, build_data.system_vectors)
end

""" Apply all four stencils to the current cell.

# Arguments
- `inz::Int`:
    - Number of non-zero matrix elements so far.

- `cell_indices::CellIndices`:
    - Cell indices of the current cell.

- `build_data::StokesBuildData`:
    - Build data.

# Returns
- `inz::Int`:
    - Updated index of non-zero matrix arrays (Lii, Ljj, Li, Lj and Lv).
"""
function apply_stencils(
    inz::Int,
    cell_indices::CellIndices,
    build_data::StokesBuildData
)::Int
    inz = StencilForVxStokesUnknown.apply_stencil(inz, cell_indices, build_data)
    inz = StencilForVyStokesUnknown.apply_stencil(inz, cell_indices, build_data)
    inz = StencilForVzStokesUnknown.apply_stencil(inz, cell_indices, build_data)
    inz = StencilForPContinuityUnknown.apply_stencil(inz, cell_indices, build_data)
    return inz
end

end # module
