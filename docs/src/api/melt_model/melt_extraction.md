# MeltExtraction

Melt is incrementally extracted from migration domains defined by local 
maxima (-y is in upward direction) along the top of the partially molten
region of the mantle. The top of the partially molten region can be 
smoothed using the parameter `smoothing_radius_drainage`. 

The total volume of melt to be extracted is calculated based on the 
volume of melt in the partially molten mantle region and the extraction 
fraction.

Extracted melt is instantaneously emplaced using one of two models:
    1. Partially Molten Emplacement (`iuse_shallow_mantle_emplacement = 0`): 
        melt is emplaced at the shallowest marker of the partially molten mantle 
        in a given migration domain.
    2. Mantle Emplacement (`iuse_shallow_mantle_emplacement = 1`): melt is 
        emplaced at the shallowest marker of the mantle located within an 
        injection width centered on the local maximum of the partially molten 
        mantle domain. The injection width is controlled by the characteristic 
        injection width parameter and will be automatically adjusted to ensure 
        emplaced magma height does not exceed a threshold (`magma_height_limit`).

To improve computational efficiency a subset of the mantle domain is
searched for the shallowest mantle marker. The mantle search width 
defines the width of this search domain, which is centered on the initial
shallowest partially molten mantle marker. The mantle search domain 
should be larger than the injection width otherwise erroneous injection 
may occur. A warning is printed if the injection domain is outside of the 
mantle search domain.

Emplacement can be controlled by a probabilistic normal distribution 
model by setting use_random_injection_subdomain to True. With this option
the injection width is divided into injection subdomains selected using 
a normal probability distribution function. The shallowest mantle marker 
within the selected subdomain is converted into magma.

The effects of fractionation in sills can be approximated by setting
`iuse_gabbroic_fractionation = 1`. This model changes the composition
of gabbroic magma to fractionated gabbroic magma once magma enters the
crustal domain. Entry into the crustal domain is controlled by the
`fractionation_threshold_limit` parameter which is the distance from the
Moho where the composition of gabbroic magma is changed to fractionated
or layered gabbroic magma. This change in composition involves an increase
in the solidus leading to regions of partially molten fractioned gabbro
as opposed to large regions of pure gabbroic magma inconsistent with
geophysical observations.

# Initialization

```@docs
MeltModel.Extraction.initialize!
```