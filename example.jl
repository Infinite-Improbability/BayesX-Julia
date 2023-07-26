include("src/BayesJ.jl")
using .BayesJ
using Unitful, DimensionfulAngles
using Profile

# ENV["JULIA_DEBUG"] = "BayesJ"

data = FITSData(
    "/home/ryan/data/chandra/4361/manual3/repro/acisf04361_repro_evt2.fits",
    "/home/ryan/data/chandra/4361/manual3/repro/bg_trimmed_300-7000.fits",
    "/home/ryan/data/chandra/4361/manual3/repro/specx/specx.arf",
    "/home/ryan/data/chandra/4361/manual3/repro/specx/specx.rmf",
    0.492u"arcsecondᵃ"
)

priors = [LogUniformPrior(1.0e14, 1.0e17), UniformPrior(0.01, 1.0)]
sample(data, (0.3:0.01:7.0)u"keV", priors, nHcol=3.89, redshift=0.164)