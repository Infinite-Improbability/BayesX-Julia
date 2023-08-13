# Setup

Using Julia's [package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/) run

```
(@v1.8) pkg> activate .
  Activating project at `~/git/BayesX-Julia`

(BayesX-Julia) pkg> instantiate
```
This project includes some dependencies not on the main Julia registry. If packages fail to install you may need to add astro-group-bristol/AstroRegistry, as explained [here](https://github.com/astro-group-bristol/AstroRegistry).
Alternatively you can locate the specific packages on GitHub and install them using the repository URLs.

After installation the sample analysis can be launched with
```sh
julia --project=@. examples/example.jl
```

# MPI Support
Ultranest supports MPI. To use launch your Julia scripts using `mpiexec` or equivalent.

While the Python dependencies will install automatically if not found, this can be cause segfaults as MPI.jl is trying to use one MPI library and mpi4py another. I suggest manually installing `ultranest` and `mpi4py` in a conda environment. This allows you to use [force Conda to use an external MPI library](https://conda-forge.org/docs/user/tipsandtricks.html#using-external-message-passing-interface-mpi-libraries).

Then use `julia --project -e 'using MPIPreferences; MPIPreferences.use_system_binary()'` to [set the MPI library used by Julia](https://juliaparallel.org/MPI.jl/stable/configuration/).

Finally you can tell PyCall to use the Python from Conda by [building it as explained here](https://docs.juliahub.com/PyCall/GkzkC/1.92.0/#Specifying-the-Python-version). In short activate the conda environment and run
```julia
julia> ENV["PYTHON"]="/path/to/python"
julia> using Pkg
julia> Pkg.precompile("PyCall")
julia> Pkg.build("PyCall")
    Building Conda ─→ `~/.julia/scratchspaces/44cfe95a-1eb2-52ea-b672-e2afdf69b78f/8c86e48c0db1564a1d49548d3515ced5d604c408/build.log`
    Building PyCall → `~/.julia/scratchspaces/44cfe95a-1eb2-52ea-b672-e2afdf69b78f/43d304ac6f0354755f1d60730ece8c499980f7ba/build.log`
Precompiling project...
  5 dependencies successfully precompiled in 10 seconds. 278 already precompiled.
```