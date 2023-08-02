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

priors = [LogUniformPrior(1.0e13, 3.0e15), UniformPrior(0.01, 1.0), UniformPrior(-10, 10), UniformPrior(-10, 10)]
sample(data, (0.7:0.01:7.0)u"keV", priors, 3.89e20u"cm^-2", 0.164)