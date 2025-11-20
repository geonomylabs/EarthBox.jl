# Marker-Grid Interpolation

Marker-in-cell methods require tensor and scalar parameters to be interpolated back and forth
between the the Eulerian grid and the Lagrangian markers. We use the first-order bilinear 
interpolation scheme for interpolating marker information to Eulerian grid nodes:

###### eq:bilinear_interp2grid
```math
S_{i,j} = 
    \frac{
        \displaystyle \sum_{m}S_m w_{m(i,j)}
    }{
        \displaystyle \sum_{m}w_{m(i,j)}
    },
\quad w_{m(i,j)} = 
    \left(
        1 - \frac{\Delta x_m}{\Delta x}
    \right)
    \left(
        1 - \frac{\Delta y_m}{\Delta y}
    \right)
```

where ``S(i, j)`` is the interpolated parameter value at the Eulerian grid node ``(i, j)``, 
``S_m`` is the parameter value at the marker ``m`` located within one of four surrounding grid 
cells, ``w_{m(i,j)}`` is the weight of the marker ``m`` at the Eulerian grid node ``(i, j)``, 
``\Delta x_m`` and ``\Delta y_m`` are the distances between the marker ``m`` and the Eulerian 
grid node ``(i, j)`` in the ``x`` and ``y`` directions, respectively, and ``\Delta x`` and 
``\Delta y`` are the grid spacings in the ``x`` and ``y`` directions, respectively. With our 
modeling approach the index i corresponds to the vertical y-axis that increases with depth and 
the index j corresponds to the x-axis. We refer to two types of interpolation based on the size 
of the search radius used in [Eq.](@ref eq:bilinear_interp2grid): (1) inclusive interpolation 
where markers within the search radius equal to ``\Delta x`` and ``\Delta y`` are used in the 
interpolation, and (2) exclusive interpolation where markers within the search radius equal to 
``0.5\Delta x`` and ``0.5\Delta y`` are used in the interpolation.

A similar first-order bilinear scheme is used for interpolating parameter values from Eulerian grid 
nodes to Lagrangian markers:

###### eq:bilinear-interp2markers
```math
\begin{split}
S_m = 
S_{i_{ul,m}, j_{ul,m}}w_{UL,m} + S_{i_{ul,m}, j_{ul,m}+1}w_{UR} 
+ S_{i_{ul,m}+1, j_{ul,m}}w_{LL} + S_{i_{ul,m}+1, j_{ul,m}+1}w_{LR}
\end{split}
```

where ``S_m`` is the interpolated parameter value at the Lagrangian marker ``m``, ``i_{ul,m}`` is 
the y-index of the upper left grid node for marker ``m``, ``j_{ul,m}`` is the x-index of the 
upper-left grid node for marker ``m``, ``S_{i_{ul,m},j_{ul,m}}``,  ``S_{i_{ul,m}, j_{ul,m}+1}``, 
``S_{i_{ul,m}+1, j}``, and ``S_{i_{ul,m}+1, j_{ul,m}+1}`` are the parameter values at the four 
surrounding Eulerian grid nodes, and ``w_{UL,m}``, ``w_{LL,m}``, ``w_{UR,m}``, and ``w_{LR,m}`` 
are the weights for the upper-left, lower-left, upper-right, and lower-right grid nodes, 
respectively. The weights of the four surrounding Eulerian grid nodes for marker ``m`` are 
calculated using the following equations:

###### eq:node_marker_weights
```math
\begin{split}
    w_{UL,m} & = (1.0 - \Delta x'_{UL})(1.0 - \Delta y'_{UL}) \\
    w_{LL,m} & = (1.0 - \Delta x'_{UL})\Delta y'_{UL}\\
    w_{UR,m} & = \Delta x'_{UL}(1.0 - \Delta y'_{UL}) \\
    w_{LR,m} & = \Delta x'_{UL}\Delta y'_{UL} \\
\end{split}
```

where ``\Delta x'_{UL}`` and ``\Delta y'_{UL}`` are the normalized distances in the ``x`` and 
``y`` directions, respectively, between the marker and the upper-left node of the cell containing 
the marker as defined by the following equations:

###### eq:normalized_upper_left_distances
```math
\begin{split}
\Delta x'_{UL} = \frac{\left(x_m - x_{b,j_{ul,m}}\right)}{\Delta x_{b,j_{ul,m}}} \\
\Delta y'_{UL} = \frac{\left(y_m - y_{b,i_{ul,m}}\right)}{\Delta y_{b,i_{ul,m}}}
\end{split}
```

where ``x_m`` and ``y_m`` are the coordinates of the marker ``m``, ``x_{b,j_{ul,m}}`` and 
``y_{b,i_{ul,m}}`` are the coordinates of the upper-left basic grid node of the cell containing 
the marker ``m``, and ``\Delta x_{b,j_{ul,m}}`` and ``\Delta y_{b,i_{ul,m}}`` are the basic grid 
spacings in the ``x`` and ``y`` directions, respectively. The upper-left indices ``i_{ul,m}`` and 
``j_{ul,m}`` of the cell containing the marker and weights ``w_{UL,m}``, ``w_{LL,m}``, 
``w_{UR,m}``, and ``w_{LR,m}`` are pre-computed for each marker using bisection to improve the 
efficiency of the interpolation process. 

These pre-computed weights and indices can also be used to optimize [Eq.](@ref eq:bilinear_interp2grid) 
by first calculating numerator and denominator terms of the bilinear average for each grid node by 
looping over each marker ``m`` and applying the following equations:

###### eq:bilinear-numerator
```math
\begin{split}
    S_{avg,n(i_{ul,m}, j_{ul,m})} & = S_{avg,n(i_{ul,m}, j_{ul,m})} + S_m w_{UL,m} \\
    S_{avg,n(i_{ul,m}, j_{ul,m}+1)} & = S_{avg,n(i_{ul,m}, j_{ul,m}+1)} + S_m w_{UR,m} \\
    S_{avg,n(i_{ul,m}+1, j_{ul,m})} & = S_{avg,n(i_{ul,m}+1, j_{ul,m})} + S_m w_{LL,m} \\
    S_{avg,n(i_{ul,m}+1, j_{ul,m}+1)} & = S_{avg,n(i_{ul,m}+1, j_{ul,m}+1)} + S_m w_{LR,m} \text{,}\\
\end{split}
```
and

###### eq:bilinear-denominator
```math
\begin{split}
    S_{avg,d(i_{ul,m}, j_{ul,m})} & = S_{avg,d(i_{ul,m}, j_{ul,m})} + w_{UL,m} \\
    S_{avg,d(i_{ul,m}, j_{ul,m}+1)} & = S_{avg,d(i_{ul,m}, j_{ul,m}+1)} + w_{UR,m} \\
    S_{avg,d(i_{ul,m}+1, j_{ul,m})} & = S_{avg,d(i_{ul,m}+1, j_{ul,m})} + w_{LL,m} \\
    S_{avg,d(i_{ul,m}+1, j_{ul,m}+1)} & = S_{avg,d(i_{ul,m}+1, j_{ul,m}+1)} + w_{LR,m} \\
\end{split}
```

where  subscript ``n`` denotes the numerator term, and subscript ``d`` denotes the denominator
term and grid arrays ``S_{avg,n}(i, j)`` and ``S_{avg,d}(i, j)`` are initialized with zero values. 
The bilinear average at each grid node ``(i, j)`` is then calculated as:

###### eq:bilinear-average
```math
S_{i, j} = \frac{S_{avg,n}(i, j)}{S_{avg,d}(i, j)} \text{.}\\
```