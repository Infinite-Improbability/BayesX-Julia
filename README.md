# BayesJ

[![Run tests](https://github.com/Infinite-Improbability/BayesX-Julia/actions/workflows/tests.yml/badge.svg)](https://github.com/Infinite-Improbability/BayesX-Julia/actions/workflows/tests.yml)
[![Documentation](https://github.com/Infinite-Improbability/BayesX-Julia/actions/workflows/documentation.yml/badge.svg)](https://github.com/Infinite-Improbability/BayesX-Julia/actions/workflows/documentation.yml)

Bayesian inference on galaxy clusters in X-ray.

Based on the original Fortran version presented in [Olamaie 2015](https://doi.org/10.1093/mnras/stu2146).

## Setup

BayesJ is an unregistered Julia package. Install it with
```julia
pkg> add https://github.com/Infinite-Improbability/BayesX-Julia
```
You will need the AstroRegistry from https://github.com/astro-group-bristol/AstroRegistry.

For details see the [documentation](https://infinite-improbability.github.io/BayesX-Julia/dev/setup/).

## Known Issues
- When starting an analysis using `mpiexec` immediately after updating SpectralFitting.jl Julia has reported an inability to find SpectralFitting with compilecache. The run has later errored between creating the new run directory and starting sampling, with the message `MethodError: Cannot convert an object of type Interpolations.Extrapolation{...} to an object of type Float64`. To resolve this issue start an analysis without MPI so the library is compiled. A call to `precompile SpectralFitting` in the Julia pkg REPL may also work.

## Testing
The units tests use the standard Julia configuration. Just run
```julia
julia> using Pkg
julia> Pkg.test("BayesJ")
```
