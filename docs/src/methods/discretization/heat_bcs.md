 # Heat Conduction Equation Boundary Conditions

Earthbox uses the approach of [gerya2010](@citet) for specifying boundary conditions  
for the heat-conduction equations on a basic grid whereby simple coefficients are used
that control type of boundary condition. 

## Heat Boundary Conditions: Top and Bottom Boundaries

Boundary conditions for the heat conduction equation are implemented along the top and bottom 
boundaries using the following equations:

###### eq:heat-bc-top-bottom
```math
\begin{split}
    T_{(1,j)_b} & = T_{top} + f_T^{top}T_{(2,j)_b} \\
    T_{(N_{y,b},j)_b} & = T_{bot} + f_T^{bot}T_{(N_{y,b}-1,j)_b} \\
\end{split}
```

where ``j`` ranges from ``1`` to ``N_{x,b}``, ``f_T^{top}`` and ``f_T^{bot}`` are the boundary 
condition coefficients for the top boundary and bottom boundaries, respectively, and ``T_{top}`` 
and ``T_{bot}`` are the prescribed temperature or gradient values at nodes along the top and 
bottom boundaries, respectively. The coefficients ``f_T^{top}`` and ``f_T^{bot}`` control how 
the temperature at boundary nodes is related to the internal temperature and can have a value 
of -1, 0 or 1 depending on the type of boundary condition.

# Heat Boundary Conditions: Left and Right Boundaries

Boundary conditions for the heat conduction equation are implemented along the left and right 
boundaries using the following equations:

###### eq:heat-bc-left-right
```math
\begin{split}
    T_{(i,1)_b} & = T_{left} + f_T^{left}T_{(i,2)_b} \\
    T_{(i,N_{x,b})_b} & = T_{right} + f_T^{right}T_{(i,N_{x,b}-1)_b} \\
\end{split}
```

where ``i`` ranges from ``1`` to ``N_{y,b}``, ``f_T^{left}`` and ``f_T^{right}`` are the 
boundary condition coefficients for the left boundary and right boundaries, respectively, and 
``T_{left}`` and ``T_{right}`` are the prescribed temperature or gradient values at nodes along 
the left and right boundaries, respectively. The coefficients ``f_T^{left}`` and ``f_T^{right}`` 
control how the temperature at boundary nodes is related to the internal temperature and can 
have a value of 0 or 1 depending on the type of boundary condition.

Boundary nodes have a large matrix coefficient of 1 for the unknown temperature ``T`` and a right 
hand side value equal to the prescribed temperature value ``T_{top}``, ``T_{bot}``, ``T_{left}`` 
or ``T_{right}`` for the top, bottom, left and right boundaries, respectively. The internal nodes 
used in equations [Eq.](@ref eq:heat-bc-top-bottom) and [Eq.](@ref eq:heat-bc-left-right) are
associated with unknowns that have coefficients equal to ``-f_T^{top}``, ``-f_T^{bot}``, 
``-f_T^{left}`` and ``-f_T^{right}`` for the top, bottom, left and right boundaries, respectively, 
after rearranging unknowns to the left hand side of equations [Eq.](@ref eq:heat-bc-top-bottom) 
and [Eq.](@ref eq:heat-bc-left-right).

# Constant Temperature Boundary Conditions

Constant temperature boundary conditions are implemented by setting ``f_T^{top}``, ``f_T^{bot}``, 
``f_T^{left}`` and ``f_T^{right}`` to 0 and setting ``T_{top}``, ``T_{bot}``, ``T_{left}`` and 
``T_{right}`` to the prescribed temperature values at the top, bottom, left and right boundaries, 
respectively. Constant gradient boundary conditions are implemented by setting ``f_T^{top}``, 
``f_T^{bot}``, ``f_T^{left}`` and ``f_T^{right}`` to 1 and setting ``T_{top}``, ``T_{bot}``, 
``T_{left}`` and ``T_{right}`` to the prescribed temperature gradient values at the top, bottom, 
left and right boundaries, respectively. 

## Zero Flux Boundary Conditions

Zero heat flux boundary conditions are implemented by setting ``f_T^{top}``, ``f_T^{bot}``, 
``f_T^{left}`` and ``f_T^{right}`` to 1 and setting ``T_{top}``, ``T_{bot}``, ``T_{left}`` and 
``T_{right}`` to zero at the top, bottom, left and right boundaries, respectively. 

## Heat Boundary Conditions for Extensional Models

For extensional models constant temperature boundary conditions are used along the top and bottom 
boundaries and zero heat flux boundary conditions are used along the left and right boundaries as 
described in the following equations:

###### eq:heat-bc-extension-models
```math
\begin{split}
    T_{(1,j)_b} & = T_{top} \\
    T_{(N_{y,b},j)_b} & = T_{bot} \\
    T_{(i,1)_b} & = T_{(i,2)_b} \\
    T_{(i,N_{x,b})_b} & = T_{(i,N_{x,b}-1)_b} \text{.} \\
\end{split}
```