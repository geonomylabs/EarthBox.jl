# Staggered Grid Generation

A grid type option referred to as T-type involves local refinement in a upper central domain
of the model and smooth variations of grid resolution in surrounding lower resolution domains. T-type
refinement is enabled by setting the `grid_type` argument equal to `:TtypeRefinedGrid` (See 
[StaggeredGrid.initialize!](@ref EarthBox.StaggeredGrid.initialize!)). 

The high-resolution domain in the T-type grid has basic grid spacings ``\Delta x_{hr}`` and 
``\Delta y_{hr}`` in the ``x-`` and ``y``-directions, respectively, and extends laterally from 
x-coordinate ``x_o`` to ``x_f`` and vertically from the top of the model domain to a depth of 
``y_f``. The x-indices of the basic grid corresponding to ``x_o`` and ``x_f`` are ``j_{x_o}`` 
and ``j_{x_f}``, respectively, and are calculated as follows:

###### eq:ttype-bounding-j-indices
```math
\begin{split}
    j_{x_o} & = int\left(\frac{N_{x,b} - N_{x,hr}}{2}\right) + 1\\
    j_{x_f} & = j_{x_o} + int\left( \frac{\left(x_f - x_o\right)}{\Delta x_{hr}} \right)
\end{split}
```

Starting with ``x_{b,j_{x_o}} = x_o`` the x-coordinates within the high resolution domain from 
``j=j_{x_o}+1`` to ``j=j_{x_f}`` are calculated as follows:

###### eq:ttype-x-coords-highres
```math
x_{b,j} = x_{b,j-1} + \Delta x_{hr} \text{.}
```

Starting with index ``j=2`` and ending with index ``j=j_{x_o}`` the x-coordinates to the left of 
``x_o`` are calculate using the following equation [gerya2010](@cite):

###### eq:ttype-x-coords-left
```math
x_{b,j} = x_{b,j-1} + \Delta x_{hr}(f_{smooth,L})^{j_{x_o} - j + 1}
```

where ``f_{smooth,L}`` is a smoothing factor that is calculated iteratively starting with an initial 
value of ``f_{smooth,L} = 1.1`` using the following equation:

###### eq:ttype-left-smoothing
```math
f_{smooth,L} = \left[
    1 
    + \frac{x_o}{\Delta x_{hr}}\left(1 - \frac{1}{f_{smooth,L}}\right) 
    \right]^{\frac{1}{N_{steps,L}}}
```

where ``N_{steps,L} = j_{x_o} - 1``. The x-coordinates to the right of ``x_f`` starting at index ``j_{x_f}+1`` 
are calculated using the following equation [gerya2010](@cite):

######  eq:ttype-x-coords-right
```math
x_{b,j} = x_{b,j-1} + \Delta x_{hr}(f_{smooth,R})^{j - j_{x_f}}
```

where ``f_{smooth,R}`` is a smoothing factor that is calculated iteratively starting with an initial
value of ``f_{smooth,R} = 1.1`` using the following equation:

###### eq:ttype-right-smoothing
```math
f_{smooth,R} = \left[
    1 
    + \frac{x_{size} - x_{b,j_{x_f}}}{\Delta x_{hr}}\left(1 - \frac{1}{f_{smooth,R}}\right) 
    \right]^{\frac{1}{N_{steps,R}}}
```

where ``N_{steps,R} = \left(N_{x,b} - j_{x_f}\right)``. 

The spacing of the basic grid in the x-direction is given by:

###### eq:x-spacing-basic
```math
\Delta x_{b,j} = x_{b,j+1} - x_{b,j} \text{.}
```

where ``j`` ranges from ``1`` to ``N_{x,b}-1``.

Starting with ``y_{b,i} = 0``, the y-coordinates of the high-resolution domain from indices 
``i=2`` to ``i=i_{y_f}`` are calculated as follows:

###### eq:ttype-y-coords-highres
```math
y_{b,i} = y_{b,i-1} + \Delta y_{hr} \text{,}
```

where ``i_{y_f}`` is the index of the basic grid corresponding to the depth ``y_f`` and is calculated as follows:

###### eq:ttype-y-index-final
```math
i_{y_f} = int\left(\frac{y_f}{\Delta y_{hr}}\right) + 1 \text{.}
```

The y-coordinates below the high-resolution domain starting at index ``i = i_{y_f}+1`` are calculated using 
the following equation [gerya2010](@cite):

###### eq:ttype-y-coords-below-hr
```math
y_{b,i} = y_{b,i-1} + \Delta y_{hr}(f_{smooth,B})^{i - i_{y_f}}
```

where ``f_{smooth,B}`` is a smoothing factor that is calculated iteratively starting with an initial
value of ``f_{smooth,B} = 1.1`` using the following equation:

###### eq:ttype-y-smoothing
```math
f_{smooth,B} = \left[
    1 
    + \frac{y_{size} - y_{b,i_{y_f}}}{\Delta y_{hr}}\left(1 - \frac{1}{f_{smooth,B}}\right) 
    \right]^{\frac{1}{\left(N_{y,b} - i_{y_f}\right)}} \text{.}
```

The spacing of the basic grid in the y-direction is given by:

###### eq:y-spacing-basic
```math
\Delta y_{b,i} = y_{b,i+1} - y_{b,i} \text{.}
```

where ``i`` ranges from ``1`` to ``N_{y,b}-1``.

The pressure grid, denoted by subscript ``p``, has nodes located at the center of basic grid cells 
with ``(N_{x,b}-1)`` nodes in the ``x``-direction and ``(N_{y,b}-1)`` nodes in the ``y``-direction 
[fig:staggered-grid](@ref fig:staggered-grid). The x- and y-coordinates of the pressure grid are ``(x_{p,j}, y_{p,i})`` 
where ``j`` is in the ``x``-direction ranging from ``1`` to ``N_{x,b}-1`` and ``i`` is in the ``y``-direction ranging 
from ``1`` to ``N_{y,b}-1``. The pressure grid spacings in the ``x``- and ``y``-directions are denoted by 
``\Delta x_{p,j}`` where ``j`` ranges from ``1`` to ``N_{x,b}-2`` and ``\Delta y_{p,i}`` where ``i`` ranges from ``1`` 
to ``N_{y,b}-2``, respectively. Pressure unknown ``P_{(i,j)_p}``, deviatoric normal stress 
``\sigma'_{xx(i,j)_p}``, deviatoric normal strain rate ``\dot \epsilon'_{xx(i,j)_p}`` and 
visco-plastic normal viscosity ``\eta_{vp(i,j)_p}`` are defined on the nodes of the pressure grid. 
The spacing of the pressure grid in the x-direction is given by:

###### eq:pressure-x-spacing
```math
\Delta x_{p,j} = 0.5\left(x_{b,j+2} - x_{b,j}\right)
```

where ``j`` ranges from ``1`` to ``N_{x,b}-2``. The spacing of the pressure grid in the y-direction is given by:

###### eq:pressure-y-spacing
```math
\Delta y_{p,i} = 0.5\left(y_{b,i+2} - y_{b,i}\right)
```

where ``i`` ranges from ``1`` to ``N_{y,b}-2``. Starting with ``x_{p,1} = 0.5(x_{b,1}+x_{b,2})``, the 
x-coordinates of the pressure grid from indices ``j=2`` to ``j=N_{x,b}-1`` are calculated as follows:

###### eq:pressure-x-coords
```math
x_{p,j} = x_{p,j-1} + \Delta x_{p,j-1} \text{.}
```

Starting with ``y_{p,1} = 0.5(y_{b,1}+y_{b,2})``, the y-coordinates of the pressure grid from indices
``i=2`` to ``i=N_{y,b}-1`` are calculated as follows:

###### eq:pressure-y-coords
```math
y_{p,i} = y_{p,i-1} + \Delta y_{p,i-1} \text{.}
```

The ``v_x`` grid, denoted by subscript ``v_x``, has nodes located at the center of left- and 
right-sides of basic-grid cells and outside the model domain in the y-direction with ``N_{x,b}`` nodes
in the ``x``-direction and ``(N_{y,b}+1)`` nodes in the ``y``-direction [fig:staggered-grid](@ref fig:staggered-grid). 
The x- and y-coordinates of the ``v_x`` grid are ``(x_{b,j}, y_{v_x,i})`` where ``j`` ranges from ``1`` to 
``N_{x,b}`` and ``i`` ranges from ``1`` to ``N_{y,b}+1``. The grid spacings in the y-direction of the ``v_x`` 
grid are denoted by ``\Delta y_{v_x,i}`` where ``i`` ranges from ``1`` to ``N_{y,b}``. Spacing in the 
x-direction is equal to ``\Delta x_{b,j}``. The unknown for the x-component of velocity ``v_{x(i,j)v_x}`` 
is defined on the nodes of the ``v_x`` grid. Starting with ``y_{v_x,1} = 0.5(y_{b,1} - 0.5\Delta y_{b,1})`` 
and ending with ``y_{v_x,N_{y,b}+1} = y_{b,N_{y,b}} + 0.5\Delta y_{b,N_{y,b}-1}``, the y-coordinates of 
the ``v_x`` grid from indices ``i=2`` to ``i=N_{y,b}`` are calculated as follows:

###### eq:vx-grid-y-coords
```math
y_{v_x,i} = 0.5(y_{b,i} + y_{b,i-1}) \text{.}
```

Starting with ``\Delta y_{v_x,1} = \Delta y_{b,1}`` and ending with ``\Delta y_{v_x,N_{y,b}} = \Delta y_{b,N_{y,b}-1}``, the grid spacings in the y-direction of the ``v_x`` grid from indices
``i=2`` to ``i=N_{y,b}`` are calculated as follows:

###### eq:vx-y-spacing
```math
\Delta y_{v_x,i} = 0.5(y_{b,i+1} - y_{b,i-1})
```

The ``v_y`` grid, denoted by subscript ``v_y``, has nodes located at the center top and bottom 
of basic-grid cells and outside the model domain in the x-direction with ``N_{x,b}+1`` nodes
in the ``x``-direction and ``(N_{y,b})`` nodes in the ``y``-direction [fig:staggered-grid](@ref fig:staggered-grid). The 
x- and y-coordinates of the ``v_y`` grid are ``(x_{j,v_y}, y_{b,i})`` where ``j`` ranges from ``1`` to 
``N_{x,b}+1`` and ``i`` ranges from ``1`` to ``N_{y,b}``. The grid spacings in the x-direction of the ``v_y`` 
grid are denoted by ``\Delta x_{j,v_y}`` where ``j`` ranges from ``1`` to ``N_{x,b}``. Spacing in the 
y-direction is equal to ``\Delta y_{b,i}``. The unknown for the y-component of velocity ``v_{y(i,j)v_y}`` 
is defined on the nodes of the ``v_y`` grid. Starting with ``x_{v_y,1} = 0.5(x_{b,1} - 0.5\Delta x_{b,1})`` 
and ending with ``x_{v_y,N_{x,b}+1} = x_{b,N_{x,b}} + 0.5\Delta x_{b,N_{x,b}-1}``, the x-coordinates 
of the ``v_y`` grid from indices ``j=1`` to ``j=N_{x,b}`` are calculated as follows:

###### eq:vy-x-coords
```math
x_{v_y,j} = 0.5(x_{b,j} + x_{b,j-1}) \text{.}
```

Starting with ``\Delta x_{v_y,1} = \Delta x_{b,1}`` and ending with 
``\Delta x_{v_y,N_{x,b}} = \Delta x_{b,N_{x,b}-1}``, the grid spacings in the x-direction of the ``v_y`` grid from indices
``j=2`` to ``j=N_{x,b}-1`` are calculated as follows:

###### eq:vy-x-spacing
```math
\Delta x_{v_y,j} = 0.5(x_{b,j+1} - x_{b,j-1})
```
