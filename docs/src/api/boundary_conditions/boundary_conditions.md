# BoundaryConditions

The `BoundaryConditions` module is used to initialize the integrated
boundary condition model for the Stokes-continuity and heat equations. An
integrated boundary condition model also manages volume conservation and marker 
recycling for inflow/outflow boundary conditions.

# Boundary Condition Models

```@eval
using EarthBox
import Markdown
options = ModelManager.BoundaryConditions.get_options()
liststr = join(["- $(options[id].option_name)" for id in keys(options)], "\n")
Markdown.parse(liststr)
```

!!! tip "quick search for boundary condition model type"
    Highlight a boundary condition model name in the list above and use `Ctl-F` + `Enter`
    to navigate to a detailed description in the **Integrated Boundary Condition Model Types** section below.

# Initialization

```@docs
ModelManager.BoundaryConditions.initialize!
```

# Example