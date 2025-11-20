# Marker Distribution

Marker coordinates ``(x_m, y_m)`` are initialized using the following equations:

###### eq:marker-coordinates
```math
\begin{split}
x_m & = j_m \Delta x_m - 0.5 \Delta x_m + (r_{j_m} - 0.5) \Delta x_m \\
y_m & = i_m \Delta y_m - 0.5 \Delta y_m + (r_{i_m} - 0.5) \Delta y_m
\end{split}
```

where ``i_m`` is the marker index in the x-direction ranging from 1 to ``N_{x,m}`` 
with ``N_{x,m}`` equal to the total number of markers in the x-direction, ``j_m`` 
is the marker index in the y-direction ranging from from 1 to ``N_{y,m}`` with 
``N_{y,m}`` equal to the total number of markers in the y-direction, 
``\Delta x_m`` is the average marker spacing in the x-direction, ``\Delta y_m`` 
is the average marker spacing in the y-direction, and ``r_{i_m}`` and ``r_{j_m}`` 
are random numbers in the range ``[0,1]``. ``N_{x,m}`` and ``N_{y,m}`` are 
calculated using average marker spacings ``\Delta x_m`` and ``\Delta y_m``and model 
domain dimensions ``x_{size}`` and ``y_{size}`` as follows:

###### eq:marker-dimensions
```math
\begin{split}
    N_{x,m} = \frac{x_{size}}{\Delta x_m} \\
    N_{y,m} = \frac{y_{size}}{\Delta y_m} \text{.}
\end{split}
```

The global marker index $m$ is related to the indices $i_m$ and $j_m$ as follows:

###### eq:global-marker-index
```math
m = (j_m - 1) N_{y,m} + i_m \text{.}
```

and the total number of markers $N_{m}$ is given by:

###### eq:markers-total-number
```math
N_{m} = N_{x,m}N_{y,m} \text{.}
```