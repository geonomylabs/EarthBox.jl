 # Heat-solver Steps

[1] Update purely conductive temperature solution ``T_{solu(i,j)_b}`` and conductive temperature 
change ``\Delta T_{(i,j)_b}`` using steps \ref{step:conductive_temperature_solution} and 
[eq.temperature-change] with temperature interpolated from markers to basic grid ``T_{(i,j)_b}`` 
as the initial condition.

[2] Update marker temperature ``T_m`` for conductive temperature change ``\Delta T_{(i,j)_b}`` and 
sub-grid thermal diffusion using [Subgrid Thermal Diffusion Steps](@ref) with ``T_{(i,j)_b}``.