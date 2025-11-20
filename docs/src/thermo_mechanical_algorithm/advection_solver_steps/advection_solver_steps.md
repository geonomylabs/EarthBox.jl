# Advection-solver Steps

[1] Update marker location ``(x_m, y_m)`` and spin ``w_m`` using a 4th order Runge-Kutta scheme
and the current velocity field ``v_{x(i,j)_{vx}}`` and ``v_{y(i,j)_{vy}}``.

[2] Update marker stress ``\sigma_{xx,m}'`` and ``\sigma_{xy,m}'`` for rotation using equation
[eq:stress-rotation](@ref).

[3] Update total marker strain ``\epsilon_m`` and marker plastic strain rate 
``\dot{\epsilon}_m^{plastic}`` and plastic strain ``\epsilon_m^{plastic}`` 
([Marker Strain Update Steps](@ref)).

[4] Correct temperature in sticky air-water markers by setting temperature equal to the 
temperature at the top of the model domain.


## Marker Strain Update Steps

[1] Update total marker strain after copying current marker strain ``\epsilon_m`` to ``\epsilon_m^o``:
```math
\epsilon_m = \epsilon_m^o + \Delta t \sqrt{\dot{\epsilon}_{xx,m}^2 + \dot{\epsilon}_{xy(m)}^2}
```

[2] Update marker plastic strain rate ``\dot \epsilon _m^{plastic}`` and plastic strain 
``\epsilon _{plastic,m}`` using [eq:plastic-strain](@ref) after copying current marker plastic strain 
to \epsilon _{plastic,m}^o``
