# Installing Custom MUMPS Libraries

Using a compiled version of MUMPS that link to system libraries can lead to
significant performance improvements when solving large systems of equations
especially for the analysis phase of the MUMPS solver and enable additional
parallel ordering methods that may not be available with the Julia binaries.

The instructions here assume (1) you are working on a Debian Linux OS (e.g. Ubuntu and 
PopOS), (2) you are working with a local copy of the EarthBox project where MUMPS has 
already been installed and (3) the EarthBox project has been put in development mode in 
the base Julia environment so EarthBox can be accessed 
via `using EarthBox` (see [EarthBox Installation for Development and Experimentation](@ref)).

If you are working with a version of EarthBox installed from the registry you can perform
a similar set of steps in the base Julia environment. This will require installing
`MUMPS.jl` in the base environment first using a similar version to the one used by EarthBox.

[1] Set the environment variable `JULIA_MUMPS_LIBRARY_PATH` to the path of the directory 
that contains the local MUMPS libraries. This should be set in your `~/.bashrc` as follows:

```bash
export JULIA_MUMPS_LIBRARY_PATH="/path/to/MUMPS_x.x.x/lib"
```

This path will be used to recompile the MUMPS package for Julia. For guidance on
compiling MUMPS and generated the MUMPS library directory see [Building MUMPS Libraries](@ref).

[2] Define the `LD_LIBRARY_PATH` path to ensure that all system libraries are found 
when `MUMPS.jl` is re-built:

```bash
export LD_LIBRARY_PATH="$JULIA_MUMPS_LIBRARY_PATH:$LD_LIBRARY_PATH"
```

[3] Source your `~/.bashrc` file:

```bash
source ~/.bashrc
```

[4] Navigate to the EarthBox project directory and make a development copy of MUMPS:

```bash
julia --project=.
```

```julia
Pkg.develop("MUMPS")
```

[5] Rebuild MUMPS.jl using the development package:

```julia
Pkg.build("MUMPS")
```

The development package looks for environment variable JULIA_MUMPS_LIBRARY_PATH.
If this variable is found it will be used to recompile MUMPS.jl using the 
system MUMPS libraries.

[6] Ensure all dependencies are resolved and instantiated:

```julia
Pkg.resolve()
Pkg.instantiate()
```

# Building MUMPS Libraries

Installing MUMPS involves downloading the MUMPS source code, copying the relevant 
Makefile.inc file from the collection for different setups, editing this file 
using desired options and library locations, and then compiling the MUMPS 
libraries using make all. However, in reality a lot of editing of the 
`Makefile.inc` and editing of `Makefiles` in the src directory is often needed to 
get compilation to work. MUMPS also uses a variety of other libraries that 
can be installed with the system package manager (see the relevant 
Makefile.inc). The AI chat assistant in Cursor is very useful for figuring out
how to do this if you allow the system to perform trial and error and iterate 
on the compilation process. See the following site to download different 
versions of the MUMPS source code including the most recent release at the 
time of writing this script (5.7.3):

- [ThirdParty-Mumps](https://coin-or-tools.github.io/ThirdParty-Mumps/)

The instructions provided here assume you are working on a Debian-based Linux
OS (e.g. Ubuntu/PopOS), MUMPS version is 5.7.3 and that core dependencies have 
been installed on the system via:

```bash
apt-get install libmetis-dev libparmetis-dev libscotch-dev libptscotch-dev libatlas-base-dev openmpi-bin libopenmpi-dev liblapack-dev libscalapack-openmpi-dev
```

After installation these libraries will be located at /usr/lib/x86_64-linux-gnu, 
and associated include files are located at /usr/include after installation.

Next the instructions in the MUMPS distribution will have to be consulted and
Makefiles will have to be edited. The required edits are very system dependent 
and will most likely require trial and error. 

However, here key edits that were made to Makefiles for MUMPS 5.7.3
in order to get compilation to work for a parallel version of MUMPS that worked
well with MUMPS.jl version 1.5.1:

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

- Special care must be taken with the so called optimized options. One 
   approach is use the edits described in [Optimization Options](@ref).

## Optimization Options

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
```