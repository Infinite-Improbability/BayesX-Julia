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
julia -O3 --project=@. src/bayes.jl
```

# MPI Support
Ultranest supports MPI. To use just launch your Julia scripts using `mpiexec` or equivalent.