# Topography Grid Update Steps

[1] Advect Eulerian topography marker chain nodes ``(x_t, y_t)`` using 4th order Runge-Kutta scheme
and the current velocity field ``v_{x(i,j)_{v_x}}`` and ``v_{y(i,j)_{v_y}}``
to obtain ``(x_t', y_t')``.

[2] Interpolate advected topography nodes ``(x_{t,adv}, y_{t,adv})`` to the
Eulerian topography grid the to obtain new y-coordinates ``y_t^{new}`` for the Eulerian grid at
x-coordinates ``x_t``.

[3] Update topography marker chain nodes ``(x_t, y_t)`` by setting ``y_t`` equal to ``y_t^{new}``.