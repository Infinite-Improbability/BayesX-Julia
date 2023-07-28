include("src/BayesJ.jl")
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

priors = [LogUniformPrior(1.0e14, 1.0e17), UniformPrior(0.01, 1.0)]
sample(data, (0.3:0.01:7.0)u"keV", priors, nHcol=3.89, redshift=0.164)