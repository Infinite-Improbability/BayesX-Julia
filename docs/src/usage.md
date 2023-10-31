# Usage

We'll going to walkthrough the construction of an analysis script like in `examples/example.jl`.

First we need to import our dependencies.
```julia
using BayesJ
using Unitful, DimensionfulAngles
```

If you want to enable debug output add
```julia
ENV["JULIA_DEBUG"] = "BayesJ"
```

Now we want to load our data. Loading data from a Chandra-style events file is recommended but you can also try [`PlaintextData`](@ref)
```julia
data = FITSData(
    "data/tng/tng_projections/tng_s67_h11_x_obs_evt.fits",
    "data/tng/tng_projections/tng_s67_h11_x_bg.fits",
    "data/tng/acisi_aimpt_cy0.arf",
    "data/tng/acisi_aimpt_cy0.rmf",
    0.492u"arcsecondᵃ"
)
```
Here we have the source and background events, the ARF and RMF and the pixel size. The fancy `ᵃ` in the units for the pixel size denotes a dimensionful angle. You can get it by typing `\^a` in the Julia REPL (or VSCode if you have the Julia extension).

Now we want to define our priors.
```julia
priors_nfw = [
    UniformPrior("x0", -100.0, 100.0),
    UniformPrior("y0", -100.0, 100.0),
    UniformPrior("MT_500", 1.0e14, 1.0e15),
    NormalPrior("fg_500", 0.13, 0.01),
    UniformPrior("c_500", 3.0, 10.0),
    DeltaPrior("a", 1.0510),
    DeltaPrior("b", 5.4905),
    DeltaPrior("c", 0.3081),
    DeltaPrior("c_500_GNFW", 1.177)
]
```
The first two priors should always be called `x0` and `y0` (or x and y). They denote the cluster centre position. Other priors can be named whatever you like. These should be one prior for every positional argument to the cluster model. Currently supported priors are [`DeltaPrior`](@ref), [`LogUniformPrior`](@ref), [`NormalPrior`](@ref),[`UniformPrior`](@ref) and [`GenericPrior`](@ref), which can use arbitary distributions.

Finally we call [`sample`](@ref).
```julia
sample(
    data, # our FITSData
    (0.3u"keV", 4.0u"keV"), # (minimum energy, maximum energy)
    Model_NFW, # model function
    priors_nfw, # list of prior objects
    0.022e22u"cm^-2", # hydrogen column density
    0.5, # redshift
    (1900, 2800), # bounds on x in physical pixels
    (1900, 2800); # bounds on y in phyiscal pixels
    # everything below is optional
    bin_size=10, # bin size of spatial pixels (along each axis)
    centre_radius=0, # binned pixel radius around centre that will be excluded
    use_interpolation=false, # use interpolation in MEKAL emission calculations. Memory hungry
    use_stepsampler=false, # use stepsampler to improve efficency
    mask="data/tng/wavedetect.reg" # path to ciao style mask file. Elli
)
```
The available models are described in [Cluster Models](@ref).