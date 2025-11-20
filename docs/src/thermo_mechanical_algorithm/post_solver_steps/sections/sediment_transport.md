# Sediment Transport Update Steps

[1] Copy the topography grid ``(x_{topo,i}, y_{topo,i})`` to ``(x_{topo,i}, y_{topo,i}^o)`` to
define the initial topography grid.

[2] Calculate the initial thickness of sediments and volcanics ``H_{sed,i}^o`` on the 
topography grid by searching for the shallowest y-coordinates ``y_{sed,shallow,i}`` and deepest 
y-coordinates ``y_{sed,deep,i}`` of sediment and volcanic markers at each x-coordinate ``x_{topo,i}`` 
of the topography grid. A low-pass filter is applied to remove high-frequency noise:

###### eq:low-pass-filter
```math
    y_i^{LP}
        = \frac{
            \displaystyle \sum_{k=0}^{2N_{smooth}}
                \left(y_{(i - N_{smooth} + k)} \right)
            }{2N_{smooth}+1}
```

where ``y_i^{LP}`` is the ``i``-th y-coordinate of either the top or bottom sediment grid grid with 
a low-pass filter applied and ``N_{smooth}`` is the number of grid points that extend to the left 
and right from the current grid point ``i`` that will be included in the running average. Initial 
sediment thickness is calculated as ``H_{sed,i}^o = y_{sed,deep,i}^{LP} - y_{sed,shallow,i}^{LP}`` .

[3] Classify the shape of the topography at each x-coordinate ``x_{topo,i}`` of the topography grid. 
The shape of the topography is classified as either a sloping segment, a flat segment, local 
minimum or local maximum.

[4] Define drainage divides ``x_{divides}`` on the topography grid at each node classified as a 
local minimum, local maximum or a flat segment.

[5] Calculate downstream distances ``D_{stream,i}`` on the topography grid as the distances 
between drainage divides located to the left and right of each topography grid node:

###### eq:downstream-distances
```math
    D_{stream,i} = x_{divides,right} - x_{divides,left}
```

[6] Calculate water depth ``W_i`` on the topography grid ``(x_{topo,i}, y_{topo,i}')`` using
the following equation:

###### eq:water-depth
```math
    W_i = y_{topo,i}^o - y_{sl} \text{.}
```

where ``y_{sl}`` is the y-coordinate of the seal level.

[7] Run [Sediment Transport Time Loop](@ref).

[8] Set y-coordinate of topography grid ``y_{topo,i}`` equal to ``y_{topo,new,i}``.


## Sediment Transport Time Loop

**For each sediment transport time step do:**

[1] Calculate transport diffusivity ``\kappa_{s,i}`` using equation 
[eq:sediment-transport-diffusivity](@ref) with water depth ``W_i`` and downstream distance
``D_{stream,i}``.

[2] Calculate the y-coordinate of topography grid after sediment transport ``y_{topo,i}^{trans}`` by
solving equation [eq:sediment-transport](@ref) with ``\kappa_{s,i}`` and ``y_{topo,i}^o`` as the 
initial condition.

[3] Calculate the shallowest y-coordinates ``y_{sticky,shallow,i}`` and deepest y-coordinates
``y_{sticky,deep,i}`` of sticky markers at each x-coordinate ``x_{topo,i}`` of the topography grid,
and apply [Eq.](@ref eq:low-pass-filter) to remove high-frequency noise.

[4] Calculate the thickness of the sticky layer ``H_{sticky,i}`` on the topography grid
using ``H_{sticky,i} = y_{sticky,shallow,i} - y_{sticky,deep,i}``.

[5] Calculate the newly deposited sediment thickness at maximum compaction state 
``\Delta H_{sed,i}^o`` using equation [eq:depo-thickness-max-compaction](@ref) with 
``y_{topo,i}^o`` and ``y_{topo,i}^{trans}``.

[6] De-compact newly deposited compacted sediment thickness ``\Delta H_{sed,i}^o`` to
obtain the newly deposited sediment thickness ``\Delta H_{sed,i}^f`` by applying the conservation
of mass and the porosity-depth relationship from equation [eq:sediment-porosity](@ref).

[7] Calculate the compacted initial sediment thickness ``H_{sed,i}^f`` at a submud-depth 
``\Delta H_{sed,i}^f``. To accomplish this, the sedimentary basin is discretized into vertical 
1D meshes, and each 1D mesh is divided into 20 cells. Average properties for equation 
[eq:sediment-porosity](@ref) including ``\phi_o``, ``\lambda_{comp}`` and ``y_{submud}^{max}`` are 
calculated for each cell using all markers contained within a given cell. For each 1D mesh
cells are compacted using the conservation of mass and the porosity-depth relationship from
equation [eq:sediment-porosity](@ref) starting with the top cell and moving downward. Cell 
compaction only occurs if the new submud depth of the bottom of the cell exceeds the average maximum 
submud depth encountered by all markers within the cell. ``H_{sed,i}^f`` is calculated
by summing the compacted thickness of all cells associated with 1D meshes at each topography
grid node.

[8] Calculate compaction displacement for cell of 1D compaction meshes by summing thickness
changes for all cells below a given node. Compaction displacement for each marker 
``\Delta y_{comp,m}`` is calculated by mapping compaction displacement from 1D meshes to 
markers assuming a linear change in displacement within the cell.

[9] Calculate the top of the y-coordinates of uncompacted sedimentary basin 
``Y_{topo,precomp,i}`` at topography grid nodes by searching for the shallowest sedimentary
markers at each x-coordinate ``x_{topo,i}`` and applying [eq:low-pass-filter](@ref) to remove 
high-frequency noise.

[10] Update y-coordinates of sediment markers ``y_m`` due to compaction displacement after 
copying ``y_m`` to ``y_m^o`` using the following equation:
\begin{equation}
    y_m = y_m^o + \Delta y_{comp,m} \quad \text{if composition of marker ``m`` is sediment}
\end{equation}

[11] Calculate the top of the y-coordinates of compacted sedimentary basin 
``y_{topo,postcomp,i}`` at topography grid nodes by searching for the shallowest sedimentary
markers at each x-coordinate ``x_{topo,i}`` and applying [Eq.](@ref eq:low-pass-filter) to remove 
high-frequency noise.

[12] Calculate the compaction displacement at the sediment surface ``\Delta y_{comp,i}``
at each node of the topography grid using the following equation:
\begin{equation}
    \Delta y_{comp,i} = y_{topo,postcomp,i} - y_{topo,precomp,i}
\end{equation}

[13] Adjust sticky air and water markers assuming a linear change in the sticky layer with
zero displacement at the top of the model domain and a maximum displacement at the top
of the sediments equal to ``\Delta y_{comp,i}``.

[14] Calculate the new y-coordinate of the topography grid ``y_{topo,new,i}`` taking into
account compaction using equation [eq:sediment-compaction-correction](@ref) with ``H_{sed,i}^o``,
``\Delta H_{sed,i}^o``, ``H_{sed,i}^f`` and ``\Delta H_{sed,i}^f``.

[15] Set ``y_{topo,i}^{o}`` equal to ``y_{topo,new,i}``.

[16] Calculate downstream distances ``D_{stream,i}`` using updated y-coordinates of 
topography grid ``y_{topo,new,i}``.

[17] Calculate water depth ``W_i`` on the new topography grid using updated y-coordinates of 
topography grid ``y_{topo,new,i}``.