# Stokes-solver Steps

[1] Use the Picard iteration method and nodal plasticity to solve for x-component of 
velocity $v_{x(i,j)_{vx}}$, y-component of velocity $v_{y(i,j)_{vy}}$ and pressure 
$P_{(i,j)_p}$ and update visco-plastic normal viscosity $\eta_{vp(i,j)_p}$, 
visco-plastic shear viscosity $\eta_{vp(i,j)_b}$, deviatoric normal stress forecast 
$\sigma_{xx(i,j)_{p2}}'$ and deviatoric shear stress forecast $\sigma_{xy(i,j)_{b2}}'$ 
(see [Picard Loop](@ref)).

[2] Interpolate visco-plastic viscosity $\eta_{vp(i,j)_b}$ to markers if marker is 
located in a cell that has undergone plastic yielding (see 
[Visco-plastic Viscosity Update for Yielding on Basic Grid](@ref)).

[3] Update model time step $\Delta t$ ensuring that displacement is less than a user 
specified fraction of the staggered grid spacing based on current velocity solutions
$v_{x(i,j)_{vx}}$ and $v_{y(i,j)_{vy}}$ and grid spacings $\Delta x_{(i,j)_{vx}}$ 
and $\Delta y_{(i,j)_{vy}}$.

[4] Update deviatoric stresses on staggered grid $\sigma_{xx(i,j)_{p2}}'$ and 
$\sigma_{xy(i,j)_{b2}}'$ using equation [eq:deviatoric-stress-forecast](@ref) with 
updated deviatoric strain rates $\dot{\epsilon}_{xx(i,j)_p}'$ and 
$\dot{\epsilon}_{xy(i,j)_b}'$, updated visco-plastic viscosity $\eta_{vp(i,j)_p}$ and 
$\eta_{vp(i,j)_b}$ and updated model time step $\Delta t$.

[5] Calculate deviatoric stress changes on staggered grid 
$\Delta \sigma_{xx(i,j)_p}'$ and $\Delta \sigma_{xy(i,j)_b}'$ using 
[eq:visco-elasto-plastic-stress-changes](@ref) with updated deviatoric stresses 
$\sigma_{xx(i,j)_{p2}}'$ and $\sigma_{xy(i,j)_{b2}}'$ and original deviatoric 
stresses $\sigma_{xx(i,j)_{p1}}'$ and $\sigma_{xy(i,j)_{b1}}'$.

[6] Update marker strain rate components $\dot{\epsilon}_{xx,m}$ and 
$\dot{\epsilon}_{xy,m}$ using updated velocity solution by interpolating 
$\dot{\epsilon}_{xx(i,j)_p}$ and $\dot{\epsilon}_{xy(i,j)_b}$ to markers using 
equation [eq:bilinear-interp2markers](@ref).

[7] Update stress-based over velocity-based strain-rate ratio for markers $R_{sr,m}$ 
using equation [eq:strain_rate_ratio](@ref) with updated grid strain rate components 
$\dot{\epsilon}_{xx,(i,j)_p}$ and $\dot{\epsilon}_{xy,(i,j)_b}$, forecasted deviatoric 
grid stresses $\sigma_{xx(i,j)_{p2}}$ and $\sigma_{xy(i,j)_{b2}}$ and grid stress 
changes $\Delta \sigma_{xx(i,j)_p}'$ and $\Delta \sigma_{xy(i,j)_b}'$.

[8] Update marker pressure $P_m$ by interpolating pressure from staggered pressure 
grid $P_{(i,j)_p}$ to markers using equation [eq:bilinear-interp2markers](@ref).

[9] Update deviatoric marker stresses $\sigma_{xx,m}'$ and $\sigma_{xy,m}'$ for stress 
changes on staggered grid $\Delta \sigma_{xx(i,j)_p}'$ and 
$\Delta \sigma_{xy(i,j)_b}'$ and sub-grid stress diffusion 
$\Delta \sigma_{xx,sg(i,j)_p}'$ and $\Delta \sigma_{xy,sg(i,j)_b}'$ (see 
[Subgrid-stress Diffusion Steps](@ref)).