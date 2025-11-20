# Example Matrix

|     | 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 | 26 |
|-----|---|---|---|---|---|---|---|---|---|---|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|----|
|     | x | y | p | x | y | p | x | y | p | x | y  | p  | x  | y  | p  | x  | y  | p  | x  | y  | p  | x  | y  | p  | x  | y  | p  |
| 0 x | C | LL | L | D |  |  |  | **UR** |  | R | LR | R |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| 1 y | UR | C | U | LR | D | D |  |  |  |  | R |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| 2 p | R | D | 1 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| 3 x | U | UL |  | C | LL | L | D |  |  |  | UR |  | R | LR | R |  |  |  |  |  |  |  |  |  |  |  |  |
| 4 y |  | U |  | UR | C | U | LR | D | D |  |  |  |  | R |  |  |  |  |  |  |  |  |  |  |  |  |  |
| 5 p |  | U |  | R | D | 0 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| 6 x |  |  |  | U | UL |  | C | **LL** | L | D |  |  |  | UR |  | R | **LR** | R |  |  |  |  |  |  |  |  |  |
| 7 y |  |  |  |  |  |  |  | 1 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| 8 p | **L** |  |  |  | U |  | R | **D** | 0 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| 9 x | L |  |  |  |  |  | **U** | **UL** |  | C | LL | L | D |  |  |  | **UR** |  | R | LR | R |  |  |  |  |  |  |
| 10 y | UL | L |  | LL |  |  |  | **U** |  | UR | C | U | LR | D | D |  |  |  |  | R |  |  |  |  |  |  |  |
| 11 p |  |  |  | L |  |  |  | **U** |  | R | D | 0 |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
| 12 x |  |  |  | L |  |  |  |  |  | U | UL |  | C | LL | L | D |  |  |  | UR |  | **R** | LR | R |  |  |  |
| 13 y |  |  |  | UL | L |  | LL |  |  |  | U |  | UR | C | U | LR | **D** | D |  |  |  |  | R |  |  |  |  |
| 14 p |  |  |  |  |  |  | L |  |  |  | U |  | R | D | 0 |  |  |  |  |  |  |  |  |  |  |  |  |
| 15 x |  |  |  |  |  |  | **L** |  |  |  |  |  | U | UL |  | C | **LL** | L | D |  |  |  | UR |  | **R** | **LR** | R |
| 16 y |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 1 |  |  |  |  |  |  |  |  |  |  |
| 17 p |  |  |  |  |  |  |  |  |  | **L** |  |  |  | U |  | R | D | 0 |  |  |  |  |  |  |  |  |  |
| 18 x |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 1 |  |  |  |  |  |  |  |  |
| 19 y |  |  |  |  |  |  |  |  |  | UL | L |  | LL |  |  |  | **U** |  | **UR** | C | U | **LR** | D | D |  |  |  |
| 20 p |  |  |  |  |  |  |  |  |  |  |  |  | L |  |  |  | **U** |  | **R** | D | 0 |  |  |  |  |  |  |
| 21 x |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 1 |  |  |  |  |  |
| 22 y |  |  |  |  |  |  |  |  |  |  |  |  | UL | L |  | LL |  |  |  | U |  | **UR** | C | U | LR | **D** | D |
| 23 p |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | L |  |  |  | U |  | **R** | D | 0 |  |  |  |
| 24 x |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 1 |  |  |
| 25 y |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | 1 |  |
| 26 p |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  | L |  |  |  | U |  | **R** | **D** | 0 |

**Table:** Sparse matrix for discretized Stokes-continuity equations where ``N_{x,b} = N_{y,b} = 4``. 
For this case ``\Delta I_{v_{xR}}`` is equal to 9. Unknowns for ``v_x``, ``v_y`` and ``P'`` are 
denoted along the horizontal and vertical matrix axes as x, y, and p, respectively. The associated 
stencil is applied to each row of the matrix where the x-Stokes stencil is applied to rows labeled 
as x, the y-Stokes stencil is applied to rows labeled as y, and the continuity stencil is applied 
to rows labeled as p. Coefficients for unknowns associated with the equation generated from a given 
stencil are applied to the appropriate column of the matrix. Locations of unknowns relative to a 
main central node denoted by C are denoted as U, D, L, R, UL, UR, LL, LR which represent up, down, 
left, right, upper-left, upper-right, lower-left, and lower-right, respectively. Coefficients with 
bold font are equal to zero since the product of the prescribed unknown and associated coefficient 
from the x-stokes-coeff equation is moved to the right-hand-side term. Central main nodes 
associated with a given stencil with coefficients with values not equal to 0 or 1 are denoted as C. 
Matrix values equal one in the central node location are associated with prescribed boundary 
unknowns. Zeros along the diagonal represent pressure coefficients that are equal to zero since 
pressure does not appear in the incompressible continuity equation. Blank spaces are zero-value 
coefficients for non-diagonal unknowns. Note that the central pressure node ``P'`` has one non-zero 
coefficient equal to 1 at matrix location ``(2,2)`` where the pressure boundary condition is 
applied.