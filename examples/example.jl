using BayesJ
using Unitful, DimensionfulAngles

ENV["JULIA_DEBUG"] = "BayesJ"

data = FITSData(
    "data/tng/tng_projections/tng_s67_h11_x_obs_evt.fits",
    "data/tng/tng_projections/tng_s67_h11_x_bg.fits",
    "data/tng/acisi_aimpt_cy0.arf",
    "data/tng/acisi_aimpt_cy0.rmf",
    0.492u"arcsecondᵃ"
)

priors_v2006 = [
    DeltaPrior("x0", 56.7),
    DeltaPrior("y0", 49.0),
    UniformPrior("n0", 0.1e-3, 40.0e-3),
    UniformPrior("n02", 0.001e-1, 6.0e-1),
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
    DeltaPrior("x0", 0),
    DeltaPrior("y0", 0),
    UniformPrior("MT_500", 1.0e14, 1.0e15),
    UniformPrior("fg_500", 0.08, 0.2),
    UniformPrior("c_500", 0, 10),
    UniformPrior("α", 0.1, 1.0),
    DeltaPrior("a", 1.0510),
    DeltaPrior("b", 5.4905),
    DeltaPrior("c", 0.3081),
    DeltaPrior("c_500_GNFW", 1.177)
]
priors_nfw = [
    UniformPrior("x0", -100.0, 100.0),
    UniformPrior("y0", -100.0, 100.0),
    # DeltaPrior("x0", 56.7),
    # DeltaPrior("y0", 49.0),
    UniformPrior("MT_500", 1.0e14, 1.0e15),
    NormalPrior("fg_500", 0.13, 0.01),
    UniformPrior("c_500", 3.0, 10.0),
    DeltaPrior("a", 1.0510),
    DeltaPrior("b", 5.4905),
    DeltaPrior("c", 0.3081),
    DeltaPrior("c_500_GNFW", 1.177)
]

sample(
    data,
    (0.3u"keV", 4.0u"keV"),
    Model_NFW,
    priors_nfw,
    0.022e22u"cm^-2",
    0.5,
    (1900, 2800),
    (1900, 2800);
    bin_size=10,
    centre_radius=0,
    use_interpolation=false,
    use_stepsampler=false,
    # mask="data/tng/wavedetect.reg",
    ultranest_run_args=(max_num_improvement_loops=3, min_num_live_points=400)
)