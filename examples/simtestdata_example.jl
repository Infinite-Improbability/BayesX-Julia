include("../src/BayesJ.jl")
using .BayesJ
using Unitful, DimensionfulAngles

# ENV["JULIA_DEBUG"] = "BayesJ"

data = PlaintextData(
    "data/simtestdata/data64by64by32.txt",
    "data/simtestdata/BG64by64by32.txt",
    "data/simtestdata/ARF_32by1.txt",
    "data/simtestdata/RMF_32by32.txt",
    (32, 64, 64),
    (32, 32),
    300e3u"s",
    300e3u"s",
    0.492u"arcsecondᵃ",
    1u"cm^2"
)

priors = [LogUniformPrior(1.0e13, 1.0e15), LogUniformPrior(0.01, 1.0)]
sample(data, range(0.3u"keV", 7u"keV", 33), priors, nHcol=3.89e20u"cm^-2", redshift=0.164, use_interpolation=false)