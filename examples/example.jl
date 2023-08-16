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

priors_v2006 = [
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
priors_eiansto = [
    # UniformPrior("x0", -10.0, 10.0),
    # UniformPrior("y0", -10.0, 10.0),
    UniformPrior("MT_200", 1.0e14, 1.0e15),
    UniformPrior("fg_200", 0.08, 0.2),
    UniformPrior("α", 0.1, 1.0),
    # DeltaPrior("a", 1.0510),
    # DeltaPrior("b", 5.4905),
    # DeltaPrior("c", 0.3081),
    # DeltaPrior("c_500_GNFW", 1.177)
]
priors_nfw = [
    UniformPrior("x0", -10.0, 10.0),
    UniformPrior("y0", -10.0, 10.0),
    UniformPrior("MT_200", 1.0e14, 1.0e15),
    UniformPrior("fg_200", 0.08, 0.2),
    UniformPrior("a", 0.3, 10.0),
    UniformPrior("b", 2.0, 15.0),
    UniformPrior("c", 0.01, 1.0),
    UniformPrior("c_500_GNFW", 0.01, 6.0)
]

sample(
    data,
    (0.7:0.01:7.0)u"keV",
    Model_Einasto,
    priors_eiansto,
    3.89e20u"cm^-2",
    0.160;
    bin_size=10,
    centre_radius=0,
    use_interpolation=false
)