# BayesJ

[![Run tests](https://github.com/Infinite-Improbability/BayesX-Julia/actions/workflows/tests.yml/badge.svg)](https://github.com/Infinite-Improbability/BayesX-Julia/actions/workflows/tests.yml)
[![Documentation](https://github.com/Infinite-Improbability/BayesX-Julia/actions/workflows/documentation.yml/badge.svg)](https://github.com/Infinite-Improbability/BayesX-Julia/actions/workflows/documentation.yml)

## Setup

Using Julia's [package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/) run

```
(@v1.8) pkg> activate .
  Activating project at `~/git/BayesX-Julia`

(BayesX-Julia) pkg> instantiate
```
This project includes some dependencies not on the main Julia registry. If packages fail to install you may need to add astro-group-bristol/AstroRegistry, as explained [here](https://github.com/astro-group-bristol/AstroRegistry).
Alternatively you can locate the specific packages on GitHub and install them using the repository URLs.

After installation the sample analysis can be launched with
```shell
julia --project=@. examples/example.jl
```

## MPI Support
Ultranest supports MPI. To use launch your Julia scripts using `mpiexec` or equivalent.

While the Python dependencies will install automatically if not found, this can be cause segfaults as MPI.jl is trying to use one MPI library and mpi4py another. I suggest manually installing `ultranest` and `mpi4py` in a conda environment. This allows you to use [force Conda to use an external MPI library](https://conda-forge.org/docs/user/tipsandtricks.html#using-external-message-passing-interface-mpi-libraries).
```shell
conda install ultranest mpi4py openmpi=x.y.*=external_*
```

Then [set the MPI library used by Julia](https://juliaparallel.org/MPI.jl/stable/configuration/).
```shell
julia --project -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
```

Finally you can tell PyCall to use the Python from Conda by [building it as explained here](https://docs.juliahub.com/PyCall/GkzkC/1.92.0/#Specifying-the-Python-version). In short activate the conda environment and run
```julia
julia> ENV["PYTHON"]="/path/to/conda/env/python"
julia> using Pkg
julia> Pkg.build("PyCall")
```

## Testing
The units tests use the standard Julia configuration. Just run
```julia
julia> using Pkg
julia> Pkg.test("BayesJ")
```