include("../src/BayesJ.jl")
using .BayesJ
using Unitful, DimensionfulAngles
using DotEnv

DotEnv.config()

# ENV["JULIA_DEBUG"] = "BayesJ"

data = FITSData(
    "$(ENV["DEBUG_DATA_PATH"])/acisf04361_repro_evt2.fits",
    "$(ENV["DEBUG_DATA_PATH"])/4361_blanksky.fits",
    "$(ENV["DEBUG_DATA_PATH"])/specx/specx.arf",
    "$(ENV["DEBUG_DATA_PATH"])/specx/specx.rmf",
    0.492u"arcsecondᵃ"
)

priors = [
    UniformPrior(-10, 10), # x
    UniformPrior(-10, 10), # y
    UniformPrior(0.1, 40.0),
    UniformPrior(0.01, 6.0),
    UniformPrior(1.0, 600.0),
    UniformPrior(1.0, 100.0),
    UniformPrior(0.1, 3.0),
    UniformPrior(0.1, 2.0),
    UniformPrior(0.1, 5.0),
    UniformPrior(0.1, 5.0), # constrain ϵ<5
    UniformPrior(100.0, 1400.0),
    UniformPrior(1.0, 25.0),
    UniformPrior(0.01, 1.0),
    UniformPrior(1.0, 250.0),
    UniformPrior(0.1, 12.0),
    UniformPrior(0.01, 5.0),
    UniformPrior(-1.0, 1.0),
    UniformPrior(0.1, 6.0),
    UniformPrior(0.1, 12.0),
]
sample(data, (0.7:0.01:2.0)u"keV", priors, 3.89e20u"cm^-2", 0.160, bin_size=10, centre_radius=4)