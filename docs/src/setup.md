# Setup

BayesJ is developed for Julia 1.9, which can be obtained from https://julialang.org/. If you're unfamiliar with Julia's package manager go read its [getting started guide](https://pkgdocs.julialang.org/v1/getting-started/) quickly.

## Basic Installation
Launch the Julia REPL with
```shell
shell> julia
```
Set the Python environment used by PyCall.
```julia
julia> ENV["PYTHON"]=""
```
Press `]` to enter the package manager.

You can (optionally) activate an environment in the current folder. This helps keep dependencies tidy across projects. 
```julia
pkg> activate .
```

Setup the required package registries with
```julia
pkg> registry add General
pkg> registry add https://github.com/astro-group-bristol/AstroRegistry
```
Install BayesJ and all dependencies* with
```julia 
pkg> add https://github.com/Infinite-Improbability/BayesX-Julia
```
 Done! You can exit the package manager with backspace or Ctrl-C

\* If you want MPI support see the details below.

## Specifying an alternate Python
We use [PyCall](https://github.com/JuliaPy/PyCall.jl) to access Ultranest. This can use your system Python or a seperate minimal distribution and is controlled by the `PYTHON` environment variable. The PyCall documentation has details [here](https://github.com/JuliaPy/PyCall.jl#specifying-the-python-version). I recommend `PYTHON=""`, which uses an internal Miniconda distribution.

## MPI Support
Ultranest supports MPI for parallelisation. To use launch your Julia scripts using `mpiexec` or equivalent.

While the Python dependencies will install automatically if not found, this can be cause segfaults as MPI.jl is trying to use one MPI library and mpi4py another. I suggest manually installing `ultranest` and `mpi4py` in a conda environment. This allows you to use [force Conda to use an external MPI library](https://conda-forge.org/docs/user/tipsandtricks.html#using-external-message-passing-interface-mpi-libraries).

### MPI Setup

Find the Conda environment used by BayesJ with `conda env list`. This works even if Julia is using an internal Miniconda installation - look for environments with a `.julia` path and no name. In the example below it is the third environment.
```shell
shell> conda env list
# conda environments:
#
base                     /trimmed/for/length
ciao                     /nfs/home/coxry/.conda/envs/ciao
                         ~/.julia/conda/3/x86_64
```
Activate the environment with
```shell
shell> conda activate name/or/path
```
Install Ultranest, mpi4py and an appropriate external MPI package, as explained in the [Conda documentation](https://conda-forge.org/docs/user/tipsandtricks.html#using-external-message-passing-interface-mpi-libraries). I use OpenMPI. `x.y` should be replaced with the approriate version number, probably `4.1` for OpenMPI.
```shell
shell> conda install ultranest mpi4py openmpi=x.y.*=external_*
```
Then [set the MPI library used by Julia](https://juliaparallel.org/MPI.jl/stable/configuration/). Run this in your shell
```shell
shell> julia --project -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
```
or this in the REPL
```julia
julia> using MPIPreferences
julia> MPIPreferences.use_system_binary()
```
## Installation for Development
If you want to edit the source then you'll want to install it a bit differently.

First clone the repository. Launch the Julia REPL inside it.

Set the Python executable (as described above)
```julia
julia> ENV["PYTHON"]=""
```
Switch to pkg mode and activate the BayesJ environment.
```julia
pkg> activate .
```
Setup the required package registries
```julia
pkg> registry add General
pkg> registry add https://github.com/astro-group-bristol/AstroRegistry
```
Finally instantiate the environment to install dependencies.
```julia
pkg> instantiate
```
If you want MPI support the same instructions apply as before.