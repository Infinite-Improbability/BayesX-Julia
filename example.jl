include("src/BayesJ.jl")
using .BayesJ, Unitful, UnitfulAngles
using Profile

# ENV["JULIA_DEBUG"] = "BayesJ"
data = FITSData(
    "/home/ryan/data/chandra/4361/manual3/repro/acisf04361_repro_evt2.fits",
    "/home/ryan/data/chandra/4361/manual3/repro/bg_trimmed_300-7000.fits",
    "/home/ryan/data/chandra/4361/manual3/repro/specx/specx.arf",
    "/home/ryan/data/chandra/4361/manual3/repro/specx/specx.rmf",
    0.492u"arcsecond"
)
priors = [UniformPrior(1.0e17, 8.0e17), UniformPrior(0.08, 0.2)]
sample(data, [0.3u"keV", 7u"keV"], priors, nHcol=3.89, redshift=0.164)