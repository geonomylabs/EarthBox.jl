
# EarthBox.jl

- [EarthBox Overview](#earthbox-overview)
- [Operating System](#operating-systems)
- [Julia Installation](#julia-installation)
- [Basic EarthBox Installation](#basic-earthbox-installation)
- [Run a Model and Plot Results](#run-a-model-and-plot-results)
- [EarthBox Installation for Development and Experimentation](#earthbox-installation-for-development-and-experimentation)
- [Configuring MPI.jl to Use System Binaries](#configuring-mpijl-to-use-system-binaries)
- [Installing Custom MUMPS Libraries](#installing-custom-mumps-libraries)

## EarthBox Overview

EarthBox is a multiphase visco-elasto-plastic marker-in-cell geodynamic modeling program that 
discretizes the Stokes-continuity and heat transport equations on a staggered grid using 
conservative finite differences and free-surface stabilization. 
Advective processes are modeled using a 4th-order Runge-Kutta scheme that includes both 
grid and sub-grid changes in temperature and deviatoric stress. A variety of geologic 
processes are integrated in the EarthBox modeling system including partial melting, 
instantaneous melt extraction and transport, gabbroic melt fractionation, melt extrusion, 
lava flow, submarine hydrothermal cooling, exothermic serpentinization and sediment 
transport with erosion, deposition, compaction and an isostatically realistic global sea level.
EarthBox uses a composite viscous flow law that includes dislocation creep, diffusion creep,
Peierls creep, the effects of viscous softening associated with grain size reduction and the 
effects of partial melting. The frictional plastic model used by EarthBox includes strain weakening
and probabilistic melt damage associated with melt transport networks above zones of melt focusing.

The goals of the EarthBox project are two-fold: (1) create an easy-to-use and easy-to-modify 
geodynamic modeling tool tuned for shared memory machines and small clusters with performance 
comparable to C++ and Fortran and (2) develop a system that can be used to build
easily extensible libraries of options and model cases. EarthBox is currently 2D with an experimental 
3D multigrid solver that will be integrated into the EarthBox system during a future release.
EarthBox is developed using the Julia programming language and includes arm's-length integration with the 
parallel direct solver MUMPS that manages the stochastic 
nature of this powerful parallel direct solver that if not properly addressed can
lead to failed simulation runs.

If you use EarthBox.jl in academic, scientific, or technical work,
please consider citing or referencing this project. While not
required by the Apache License 2.0, citations help support and
acknowledge ongoing development.


Suggested citation:

Kneller, E. A. (2025). *EarthBox.jl: A visco-elasto-plastic geodynamic
modeling framework in Julia with marine and terrestrial sediment
transport, melt intrusion and extrusion, and lava flow modeling.*
https://github.com/eakneller/EarthBox.jl


## Operating Systems

Although the core functionality of Julia and EarthBox will work on Windows
and MacOS operating systems, we recommend using a Linux operating system since 
Linux yields better performance for scientific software, and the installation 
of complex scientific software often involves less overhead.

EarthBox was developed and tested on Debian-based Linux distributions like 
Ubuntu and PopOS. Therefore, you will see a strong bias toward Debian-based package 
management and the Bash shell in the EarthBox documentation. We recommend using Ubuntu
or PopOS since EarthBox was extensively tested on these operating systems. 

If you are working on Windows, we recommend using VirtualBox to build a virtual Linux
machine. See [VirtualBox](https://www.virtualbox.org).


## Julia Installation

The easiest way to get started with Julia is to use `juliaup` which involves
executing the following command on Linux systems:

```bash
curl -fsSL https://install.julialang.org | sh
```

juliaup will automatically customize your `~/.bashrc` file for using Julia.

See [Installing Julia](https://julialang.org/install/) for more details. 

Alternatively, you can install Julia directly from the latest tarball. For example:

```bash
wget https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.1-linux-x86_64.tar.gz
tar zxvf julia-1.12.1-linux-x86_64.tar.gz
```

See [Julia Tarball](https://julialang.org/downloads/platform/#linux_and_freebsd) for
more details. Custom installations of Julia will require adding the path to the Julia 
bin in your `~/.bashrc` file:

``` bash
export PATH="path/to/the/unpacked/julia/bin:$PATH"
```

and then sourcing your `~/.bashrc`:

```bash
source ~/.bashrc
```

It is often useful to define a symbolic link to the Julia directory, which includes 
the version number in the directory name, to eliminate the need to update your 
`~/.bashrc` when a new version is installed. For example,

```bash
ln -sfn julia-1.11.2 julia
```

You can customize your Julia setup by defining the following environment variables
in your `~/.bashrc`:

```bash
# Tell Julia to use all available threads
export JULIA_NUM_THREADS=auto
# Define where packages, artifacts etc. will be stored
export JULIA_DEPOT_PATH="/path/to/your/custom/julia-depot"
export JULIA_SCRATCH="/path/to/your/custom/julia-scratch"                   
export JULIA_ARTIFACTS_DIR="/path/to/your/custom/julia-artifacts"
```

When using MUMPS it is often a good idea to restrict the number of OMP threads
to 1 in your `~/.bashrc` especially if you are using a locally compiled version:

```bash
export OMP_NUM_THREADS=1
```


## Basic EarthBox Installation

You can install the registered version of EarthBox via the Julia REPL:

```julia
julia> ]
(@v1.x) pkg> add EarthBox
```

## Run a Model and Plot Results

Download the `models` directory from the EarthBox package, navigate
to `models/lithospheric_extension/extension_weak_fault`, change the 
`ROOT_OUTPUT_PATH` in `Model.jl` to a desired output location and then 
run the model using the following command:

```bash
julia Model.jl case_name=case1
```

Alternatively, you can use the REPL to run the model:

```julia
julia> include("Model.jl")
julia> Model.run_case(case_name="case1")
```

Plot model marker output files using:

```bash
julia Plot.jl marker_plots case_name=case1 istart=1 iend=10
```

Alternatively, you can plot marker output files using the REPL:

```julia
include("Plot.jl")
Plot.marker_plots(
    model_output_path = "/path/to/model_output",
    istart = 1,
    iend = 10
);
```

## Installing MPI in Base Julia Environment

You will need to install MPI in your base Julia environment, configure MPI using 
MPIPreferences and install `mpiexecjl` if you plan on using the MUMPS solver with 
EarthBox (see [Stable MPI Configuration](https://juliaparallel.org/MPI.jl/stable/usage/)):

1) Add MPI to your base environment:
    ```julia
    julia> ]
    # Use the same version of MPI used by EarthBox
    (@v1.12) pkg> add MPI@0.20.22
    ```

2) Add MPIPreferences to your base environment:
    ```julia
    julia> ]
    # Use the same version of MPIPreferences used by EarthBox
    (@v1.12) pkg> add MPIPreferences@0.1.11
    ```

3) Configure MPI to use the Julia OpenMPI binary:
    ```julia
    julia> using MPIPreferences
    julia> MPIPreferences.use_jll_binary("OpenMPI_jll")
    ```
    *An error message will be generated instructing you to restart Julia. This is OK, the 
    task was completed properly. Just restart the Julia REPL.*

4) Install `mpiexecjl` using the following:

    ```julia
    julia> using MPI
    julia> MPI.install_mpiexecjl()
    ```
    *Make a note of the path to the directory containing `mpiexecjl` produced as screen 
    output since you will add this path to your `.bashrc`*

5) Add the path to `mpiexecjl` to your `~/.bashrc`

    ```bash
    export PATH="/path/to/mpiexecjl/bin:$PATH"
    ```
    ```bash
    source ~/.bashrc
    ```


## EarthBox Installation for Development and Experimentation

The following instructions are for installing the EarthBox project for local 
experimentation that involves modifying the code:

1) Download the EarthBox project directory. 

2) Navigate to the EarthBox directory and ensure that `Manifest.toml` and 
`LocalPreferences.toml` files are deleted. Then enter the REPL with EarthBox as 
the active project:
    
    ```bash
    julia --project=.
    ```

3) From inside the REPL instantiate the EarthBox project:
    
    ```
    julia> using Pkg
    julia> Pkg.instantiate()
    ```
    This may take some time to complete especially if you are starting with a fresh
    depot.

4) With the active EarthBox project, configure MPI to use Julia binaries for OpenMPI.
    
    ```julia
    julia> using MPIPreferences
    julia> MPIPreferences.use_jll_binary("OpenMPI_jll")
    ```

    *MPI.jl should also be installed in the base Julia environment (see 
    [Installing MPI in Base Julia Environment](#installing-mpi-in-base-julia-environment)). 
    These MPI configuration steps are required for the use of MUMPS.jl with EarthBox.*

5) Add EarthBox to your current default Julia environment for development.
First exit out of the REPL if you are in it (you can use `exit()` from the REPL). 
Then re-enter the REPL with the base Julia environment and add the EarthBox 
project as a development project:

    ```julia
    julia> ]
    (@v1.x) pkg> dev /path/to/EarthBox
    ```

    This will allow you to access EarthBox via

    ``` julia
    using EarthBox
    ```

After these steps are completed you will have a basic working version of EarthBox
that uses Julia binaries for MPI and MUMPS.

Parallel analysis with PT-Scotch and ParMetis for the MUMPS solver may not 
work on your system when using Julia binaries. You may need to set the MUMPS 
`analysis_method` to SERIAL to work around these issues. However, having parallel 
analysis can significantly increase MUMPS performance especially with large 
systems of equations so it is worth the effort to get these features working. 
See [Configuring MPI.jl to Use System Binaries](#configuring-mpijl-to-use-system-binaries) 
and [Installing Custom MUMPS Libraries](#installing-custom-mumps-libraries) 
for guidance on using system MPI and a locally compiled 
version of MUMPS on Debian-type systems (e.g. Ubuntu, PopOS).


## Configuring MPI.jl to Use System Binaries

Performance can often be improved by configuring MPI.jl to use system MPI 
libraries. This section provides guidance on how to do this on Debian Linux 
operating systems (e.g. Ubuntu/PopOS) assuming that `MPI.jl` and `MPIPreferences.jl` 
have already been installed for EarthBox and the base Julia environment: 

1) Install required System libraries (OpenMPI)

    A system-level OpenMPI MPI implementation can be installed on Debian Linux systems 
    (Ubuntu/PopOS systems) using the following command:

    ```bash
    sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev libopenmpi3
    ```

2) Use `MPIPreferences` to configure `MPI.jl` to use the system binary:

    If you are working with a development version of the EarthBox project, 
    navigate to the project directory and run the REPL with the project activated:

    ```bash
    julia --project=.
    ```

    Now change the MPI configuration to use system binaries:

    ```julia
    using MPIPreferences
    MPIPreferences.use_system_binary()
    ```

    Which will generate output that looks like the following:

    ```julia
    ┌ Info: MPI implementation identified
    │   libmpi = "libmpi"
    │   version_string = "Open MPI v4.1.2, package: Debian OpenMPI, ident: 4.1.2, repo rev: v4.1.2, Nov 24, 2021\0"
    │   impl = "OpenMPI"
    │   version = v"4.1.2"
    └   abi = "OpenMPI"
    ┌ Info: MPIPreferences changed
    │   binary = "system"
    │   libmpi = "libmpi"
    │   abi = "OpenMPI"
    │   mpiexec = "mpiexec"
    │   preloads = Any[]
    └   preloads_env_switch = nothing
    ```

3) Re-instantiate your project:

    ```julia
    Pkg.resolve()
    Pkg.instantiate()
    ```

## Installing Custom MUMPS Libraries

Using a compiled version of MUMPS that links to system libraries can lead to
significant performance improvements when solving large systems of equations
especially for the analysis phase of the MUMPS solver and enable additional
parallel ordering methods that may not be available with the Julia binaries.

The instructions here assume (1) you are working on a Debian Linux OS (e.g. Ubuntu and 
PopOS), (2) you are working with a local copy of the EarthBox project where MUMPS has 
already been installed and (3) the EarthBox project has been put in development mode in 
the base Julia environment so EarthBox can be accessed 
via `using EarthBox` (see 
[EarthBox Installation for Development and Experimentation](#earthbox-installation-for-development-and-experimentation)).

If you are working with a version of EarthBox installed from the registry you can perform
a similar set of steps in the base Julia environment. This will require installing
`MUMPS.jl` in the base environment first using a similar version to the one used by EarthBox.

Use the following steps to configure `MUMPS.jl` to use MUMPS libraries that have been compiled
on your machine:

1) Set the environment variable `JULIA_MUMPS_LIBRARY_PATH` to the path of the directory 
that contains the local MUMPS libraries. This should be set in your `~/.bashrc` as follows:

    ```bash
    export JULIA_MUMPS_LIBRARY_PATH="/path/to/MUMPS_x.x.x/lib"
    ```

    This path will be used to recompile the MUMPS package for Julia. For guidance on
    compiling MUMPS and generated the MUMPS library directory see [Building MUMPS Libraries](#building-mumps-libraries).

2) Define the `LD_LIBRARY_PATH` path to ensure that all system libraries are found 
when `MUMPS.jl` is re-built:

    ```bash
    export LD_LIBRARY_PATH="$JULIA_MUMPS_LIBRARY_PATH:$LD_LIBRARY_PATH"
    ```

3) Source your `~/.bashrc` file:

    ```bash
    source ~/.bashrc
    ```

4) Navigate to the EarthBox project directory and make a development copy of MUMPS:

    ```bash
    julia --project=.
    ```

    ```julia
    Pkg.develop("MUMPS")
    ```

5) Rebuild MUMPS.jl using the development package:

    ```julia
    Pkg.build("MUMPS")
    ```

    The development package looks for environment variable JULIA_MUMPS_LIBRARY_PATH.
    If this variable is found it will be used to recompile MUMPS.jl using the 
    system MUMPS libraries.

6) Ensure all dependencies are resolved and instantiated:

    ```julia
    Pkg.resolve()
    Pkg.instantiate()
    ```

### Building MUMPS Libraries

Installing MUMPS involves downloading the MUMPS source code, copying the relevant 
Makefile.inc file from the collection for different setups, editing this file 
using desired options and library locations, and then compiling the MUMPS 
libraries using make all. However, in reality a lot of editing of the 
`Makefile.inc` and editing of `Makefiles` in the src directory is often needed to 
get compilation to work. MUMPS also uses a variety of other libraries that 
can be installed with the system package manager (see the relevant 
Makefile.inc). The AI chat assistant in Cursor is very useful for figuring out
how to do this if you allow the system to perform trial and error and iterate 
on the compilation process. See [ThirdParty-Mumps](https://coin-or-tools.github.io/ThirdParty-Mumps/) 
to download different versions of the MUMPS source code including the most recent 
release at the time of writing this document (5.7.3).

The instructions provided here assume you are working on a Debian-based Linux
OS (e.g. Ubuntu/PopOS) with `MUMPS 5.7.3` and that core dependencies have 
been installed on the system via:

```bash
apt-get install libmetis-dev libparmetis-dev libscotch-dev libptscotch-dev libatlas-base-dev openmpi-bin libopenmpi-dev liblapack-dev libscalapack-openmpi-dev
```

After installation these libraries will be located at `/usr/lib/x86_64-linux-gnu`, 
and associated include files are located at `/usr/include`.

Next the instructions in the MUMPS distribution will have to be consulted and
`Makefile.inc` will have to be configured for your system. The required 
configurations are very system dependent and will most likely require trial and 
error. 

However, here are key edits that were made to configuration files for `MUMPS 5.7.3`
in order to get compilation to work for a parallel version of MUMPS that worked
well with `MUMPS.jl` version `1.5.1`:

- Uncomment the parallel library option setup in `Makefile.inc` and 
   comment out the serial library option setup if parallel capability is 
   desired.

- Set `LIBEXT` to `.so` in the `Makefile.inc` file for shared instead of static 
   libraries which is required for the `MUMPS.jl` package.

- Set `PATH_OPT = -Wl,-rpath,$(topdir)/lib` in `Makefile.inc`. This option embeds
   the runtime library search path directly into the executable or shared 
   library being built. This means the system will know where to find the 
   required libraries at runtime without needing to modify the `LD_LIBRARY_PATH` 
   environment variable.

- Special care must be taken with the so-called optimized options. One 
   approach is to use the edits described in [Optimization Options](#optimization-options).

### Optimization Options

The following optimization options in the `Makefile.inc` led to successful compilation:

```makefile
OPTF    = -O -fopenmp -fallow-argument-mismatch
OPTL    = -O -fopenmp
OPTC    = -O -fopenmp
```

However, this required also editing the `Makefile` in the src directory by changing
`$(FPIC)` to `$(FPIC_OPT)` in the following lines starting at line 450:

```makefile
.SUFFIXES: .c .F .o
.F.o:
$(FC) $(OPTF) $(FPIC_OPT) -I. -I../include $(INCS) $(IORDERINGSF) $(ORDERINGSF) -c $*.F $(OUTF)$*.o
.c.o:
$(CC) $(OPTC) $(FPIC_OPT) -I../include $(INCS) $(CDEFS) $(IORDERINGSC) $(ORDERINGSC) -c $*.c $(OUTC)$*.o

$(ARITH)mumps_c.o:	mumps_c.c
$(CC) $(OPTC) $(FPIC_OPT) -I../include $(INCS) $(CDEFS) -DMUMPS_ARITH=MUMPS_ARITH_$(ARITH) \
        $(IORDERINGSC) $(ORDERINGSC) -c mumps_c.c $(OUTC)$@
```

Alternatively, this could have been accomplished in the `Makefile.inc` by setting 
the following:

```makefile
OPTF    = -O -fopenmp -fallow-argument-mismatch $(FPIC_OPT)
OPTC    = -O -fopenmp $(FPIC_OPT)