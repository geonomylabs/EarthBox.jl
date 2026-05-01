# Interactive tests

Standalone scripts that exercise individual EarthBox components — meant for
hands-on inspection during development, not for the automated `test/runtests.jl`
suite. Each file defines a `module` with a `run_test()` entry point and prints
results (or writes plot artifacts alongside the script).

## Run one from the shell

```bash
julia --project=. examples/interactive_tests/serpentinization/SerpentinizationTest.jl
```

The trailer

```julia
if abspath(PROGRAM_FILE) == @__FILE__
    SerpentinizationTest.run_test()
end
```

at the bottom of every demo fires only when the file is invoked directly, so
running the file as a script just works.

## Run one from the REPL

```julia
using EarthBox
include("examples/interactive_tests/serpentinization/SerpentinizationTest.jl")
SerpentinizationTest.run_test()
```

The trailer is dormant under `include`, so `run_test()` is called explicitly.
Re-`include` the same file to pick up edits.

## Layout

Subdirectories group demos by component (`compaction/`, `gravity/`,
`melt_model/`, `surface_processes/`, ...). Each leaf `*.jl` is a self-contained
demo whose module name matches the filename. PNGs sitting next to a demo are
prior outputs left in for reference.

## Adding a new demo

1. Place `MyThingTest.jl` under the appropriate component subdirectory.
2. Wrap the body in `module MyThingTest ... end`, define `run_test()::Nothing`,
   and import what you need from the public `EarthBox.*` API.
3. Append the script-mode trailer:

   ```julia
   if abspath(PROGRAM_FILE) == @__FILE__
       MyThingTest.run_test()
   end
   ```

No registration step — the directory listing is the index.
