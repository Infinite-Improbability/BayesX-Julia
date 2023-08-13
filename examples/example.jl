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
    UniformPrior("x0", -10, 10), # x
    UniformPrior("y0", -10, 10), # y
    UniformPrior("n0", 0.1, 40.0),
    UniformPrior("n02", 0.01, 6.0),
    UniformPrior("rc", 1.0, 600.0),
    UniformPrior("rc2", 1.0, 100.0),
    UniformPrior("α", 0.1, 3.0),
    UniformPrior("β", 0.1, 2.0),
    UniformPrior("β2", 0.1, 5.0),
    UniformPrior("ϵ", 0.1, 5.0), # constrain ϵ<5
    UniformPrior("rs", 100.0, 1400.0),
    UniformPrior("T0", 1.0, 25.0),
    UniformPrior("Tmin/T0", 0.01, 1.0),
    UniformPrior("rcool", 1.0, 250.0),
    UniformPrior("acool", 0.1, 12.0),
    UniformPrior("rt", 0.01, 5.0),
    UniformPrior("a", -1.0, 1.0),
    UniformPrior("b", 0.1, 6.0),
    UniformPrior("c", 0.1, 12.0),
]
sample(
    data,
    (0.7:0.01:2.0)u"keV",
    Model_Vikhlinin2006,
    priors,
    3.89e20u"cm^-2",
    0.160;
    bin_size=10,
    centre_radius=4
)