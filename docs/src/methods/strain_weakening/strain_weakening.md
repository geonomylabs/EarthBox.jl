# Strain Weakening and Melt Damage

Experimental works shows that the yield stress in equation [eq:yield-stress](@ref) decreases as 
phyllosilicate minerals and foliation form and the pressure-solution deformation mechanism becomes 
dominant as grain size is reduced in shear zones [bos02](@citep). This weakening of frictional-plastic
strength is modeled by reducing the cohesion and friction angle with increasing plastic strain as 
described by the following equations:

###### eq:strain-weakening
```math
\begin{split}
    \sigma_{c,m} = 
        \sigma_{c,m}^o
        + \frac{(\sigma_{c,m}^f - \sigma_{c,m}^o)}{(\epsilon^f - \epsilon^o)}(\epsilon_{plastic,m} - \epsilon^o)  \\
    \theta_m = 
        \theta_m^o
        + \frac{(\theta_m^f - \theta_m^o)}{(\epsilon^f - \epsilon^o)}(\epsilon_{plastic,m} - \epsilon^o)
\end{split}
```

where ``m`` is the index of the marker, ``\epsilon_{plastic,m}`` is plastic strain of the marker, 
``\epsilon^o`` and ``\epsilon^f`` are initial and final reference plastic strain, ``\sigma_{c,m}^o`` and 
``\sigma_{c,m}^f`` are initial and final reference cohesion, and ``\theta_m^f`` is the final reference 
friction angle and ``\theta_m^o`` is the final initial reference friction angle that is randomized at 
the beginning of each time step using the following equation from [naliboff17](@cite):

###### eq:friction-angle-randomization
```math
\theta_m^o = \theta'_m{^o} + (0.5 - r)f_{random}
```

where ``\theta'_{m}{^o}`` is the initial reference friction angle of the material associated 
with marker ``m``, ``r`` is a random number in the range [0,1], and ``f_{random}`` is a randomization 
factor controlling the magnitude of the perturbation and set equal to 10 for the experiments of this work.
The randomization of friction angles described in equation [Eq.](@ref eq:friction-angle-randomization) is 
used to model the evolving variability of frictional-plastic strength due to fluid-rock interactions 
[naliboff17](@citep). The final reference friction angle ``\theta_m^f`` is defined as:

###### eq:final-friction-angle
```math
\theta_m^f = \frac{\theta'_m{^f}}{\theta'_{m}{^o}}\theta_{m}^o
```

where ``\theta'_m{^f}`` is the unmodified final reference friction angle of the material associated 
with marker ``m``. Marker plastic strain ``\epsilon _{plastic,m}`` in equation [Eq.](@ref eq:strain-weakening) 
is calculate using the following equations:

###### eq:plastic-strain
```math
\begin{split}
    \sigma_{II,m}' = \sqrt{
        \left(\sigma_{xx,m}'\right)^2 + \left(\sigma_{xy,m}'\right)^2}
        \text{,} \\
    \dot \epsilon _{plastic,m} = \frac{1}{2}\sigma_{II,m}'
        \left( \frac{1}{\eta_{vp,m}} - \frac{1}{\eta_{creep,m}}\right)
        \text{,} \\
    \epsilon _{plastic,m} 
        = \epsilon _{plastic,m}^o 
            + \left( \dot \epsilon _{plastic,m} \Delta t - \dot \epsilon_{healing} \Delta t \right)             
\end{split}
```

where ``\epsilon _{plastic,m}^o`` is the current marker plastic strain, ``\dot \epsilon _{plastic,m}``
is the plastic strain rate, ``\sigma_{II,m}'`` is the second invariant of the stress tensor,
and ``\dot \epsilon_{healing}`` is the plastic healing rate.

As melt is transported through the lithosphere and crust it can weaken the bulk frictional-plastic 
yield strength of the rock by forming reactive channelized networks, increasing fluid pressure and 
forming dikes and sills that may be weaker than the surrounding rock prior to solidification. In this 
work, weakening associated with melt transport is assumed to be probabilistically localized above local 
maxima of the partially molten domain in the mantle that act as focusing points for melt migration. 
The maximum probability of melt damage is assumed to scale with amount of melt being generated in a 
melt drainage basin ``k``, which is defined by the characteristic magmatic crust height ``H_{mc,k}`` 
calculated during melt extraction (see equation [eq:char-magmatic-crust-height](@ref)). The weakening from 
damage associated with melt transport is modeled by reducing cohesion and friction angle based on a 
melt damage factor ``f_{damage,m}`` as described by the following equations:

###### eq:melt-damage-weakening
```math
\begin{split}
    \sigma_{c,m} = \frac{\sigma_{c,m}}{f_{damage,m}} \text{,}\\
    \theta_m = \frac{\theta_m}{f_{damage,m}} \text{.}
\end{split}
```

The melt damage factor ``f_{damage,m}`` is defined in a binary manner where ``f_{damage,m}`` is set equal
to 1 outside of the melt damage zone and set equal to either 1 or a constant value greater than 1 based
on a probability distribution within the melt damage zone. A melt damage zone is defined for a given
melt drainage basin ``k`` by the x-coordinate limits ``x_{damage}^{min}`` and ``x_{damage}^{max}``
and y-coordinate limit ``y_{damage}^{max}`` defined as follows:

###### eq:melt-damage-zone-limits
```math
\begin{split}
    x_{damage}^{min} & = x'_{shallow,pm,k} - \frac{W_{damage}}{2} \text{,}\\
    x_{damage}^{max} & = x'_{shallow,pm,k} + \frac{W_{damage}}{2} \text{,}\\
    y_{damage}^{max} & = y'_{shallow,pm,k}\\
\end{split}
```

where ``x'_{shallow,pm,k}`` is the average x-coordinate of the shallowest partially
molten mantle marker and ``y'_{shallow,pm,k}`` is the average y-coordinate of the shallowest partially
molten mantle marker calculated during melt extraction, and ``W_{damage}`` is the width of the
damage zone. The melt model is applied to each marker ``m`` in the melt damage zone excluding sticky
air and water markers.

The melt damage factor ``f_{damage,m}`` for a given marker ``m`` is calculated using the following
equation:

###### eq:melt-damage-factor
```math
f_{damage,m} = 
\begin{cases}
    f_{damage,max} & \text{if } r < p_{damage,m} \\
    1.0 & \text{if } r \geq p_{damage,m}
\end{cases}
```

where ``f_{damage,max}`` is the maximum melt damage factor set equal to 10 for the experiments of 
this work, and ``r`` is a random number between 0 and 1. The damage probability ``p_{damage,m}`` is 
calculated using the following equation:

###### eq:melt-damage-probability
```math
p_{damage,m} = 
\begin{cases}
    0.0 & \text{if } x_{damage}^{max} < x_m < x_{damage}^{min} \\
    \left[\frac{
            p_{damage}^{central}
        }{
            2\cos\left(\frac{\left(x_m - x'_{shallow,pm,k}\right)}{W_{damage}}2\pi\right) + 1
        } \right]
        & \text{if } x_{damage}^{min} \leq x_m \leq x_{damage}^{max}
\end{cases}
```

where the central damage zone damage probability ``p_{damage}^{central}`` is given by:

###### eq:central-damage-probability
```math
p_{damage}^{central} = 
\begin{cases}
    0.0 & \text{if } H_{mc,k} \leq H_{mc,min} \\

    p_{damage}^{max} & \text{if } H_{mc,k} \geq H_{mc,max} \\

    \frac{p_{damage}^{inter}}{H_{mc,inter} - H_{mc,min}}\Delta H_1
            & \text{if } H_{mc,min} < H_{mc,k} \leq H_{mc,inter} \\

    p_{damage}^{inter}
        + \frac{p_{damage}^{max} - p_{damage}^{inter}}{H_{mc,max} - H_{mc,inter}}\Delta H_2
                & \text{if } H_{mc,inter} < H_{mc,k} < H_{mc,max}
\end{cases}
```

where ``\Delta H_1 = (H_{mc,k} - H_{mc,min})``, ``\Delta H_2 = (H_{mc,k} - H_{mc,inter})``,
``H_{mc,k}`` is characteristic magmatic crust height for drainage basin ``k`` calculated 
during melt extraction, ``H_{mc,min}`` is the minimum reference characteristic magmatic crust height,
``H_{mc,max}`` is the maximum reference  characteristic magmatic crust height, ``H_{mc,inter}`` is the
intermediate reference characteristic magmatic crust height, ``p_{damage}^{max}`` is the maximum 
damage probability, and ``p_{damage}^{inter}`` is the intermediate damage probability.