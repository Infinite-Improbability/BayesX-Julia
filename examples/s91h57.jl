using BayesJ
using Unitful, DimensionfulAngles

data = FITSData(
    "data/tng/tng_s67_h11_y_highEres_obs_evt.fits",
    "data/tng/tng_s67_h11_y_highEres_bg.fits",
    "data/tng/acisi_aimpt_cy0.arf",
    "data/tng/acisi_aimpt_cy0.rmf",
    0.492u"arcsecondᵃ"
)

priors_v2006 = [
    NormalPrior("x0", 0.0, 20.0),
    NormalPrior("y0", 0.0, 20.0),
    LogUniformPrior("n0", 0.1e-3, 40.0e-3),
    LogUniformPrior("n02", 0.001e-2, 6.0e-1),
    LogUniformPrior("rc", 1.0, 600.0),
    LogUniformPrior("rc2", 1.0, 100.0),
    LogUniformPrior("α", 0.01, 3.0),
    LogUniformPrior("β", 0.1, 2.0),
    LogUniformPrior("β2", 0.1, 5.0),
    LogUniformPrior("ϵ", 0.1, 5.0), # constrain ϵ<5
    LogUniformPrior("rs", 100.0, 1400.0),
    LogUniformPrior("T0", 1.0, 50.0),
    LogUniformPrior("Tmin/T0", 0.01, 2.0),
    LogUniformPrior("rcool", 1.0, 250.0),
    LogUniformPrior("acool", 0.1, 12.0),
    LogUniformPrior("rt", 0.01, 10.0),
    UniformPrior("a", -1.0, 1.0),
    LogUniformPrior("b", 0.1, 6.0),
    LogUniformPrior("c", 0.1, 12.0),
]

sample(
    data,
    (0.3u"keV", 7.0u"keV"),
    Model_Vikhlinin2006,
    priors_v2006,
    0.022e22u"cm^-2",
    0.1,
    (1400, 3500),
    (1400, 3500);
    spatial_bin_size=10,
    energy_bin_size=5,
    centre_radius=0,
    use_interpolation=false,
    use_stepsampler=false,
)

