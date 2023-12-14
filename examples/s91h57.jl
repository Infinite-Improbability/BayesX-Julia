using BayesJ
using Unitful, DimensionfulAngles

data = FITSData(
    "data/tng/tng_s67_h11_y_highEres_obs_evt.fits",
    "data/tng/tng_s67_h11_y_highEres_bg.fits",
    "data/tng/acisi_aimpt_cy0.arf",
    "data/tng/acisi_aimpt_cy0.rmf",
    0.492u"arcsecondᵃ"
)

priors_nfw = [
    # UniformPrior("x0", -100.0, 100.0),
    # UniformPrior("y0", -100.0, 100.0),
    # center positions were obtained by setting other priors to reasonable fixed values then fitting center priors
    DeltaPrior("x0", 11.23),
    DeltaPrior("y0", 19.48),
    UniformPrior("MT_200", 1.0e14, 3.0e14),
    NormalPrior("fg_200", 0.108, 0.001),
    DeltaPrior("c_200", 6.0774),
    DeltaPrior("a", 1.0510),
    DeltaPrior("b", 5.4905),
    DeltaPrior("c", 0.3081),
    DeltaPrior("c_500_GNFW", 1.177)
]

priors_v2006 = [
    DeltaPrior("x0", 11.23),
    DeltaPrior("y0", 19.48),
    LogUniformPrior("n0", 0.1e-3, 40.0e-3),
    LogUniformPrior("n02", 0.001e-2, 6.0e-1),
    DependentLogUniformPrior("rc", "rc2", 600.0),
    LogUniformPrior("rc2", 1.0, 100.0),
    LogUniformPrior("α", 0.01, 3.0),
    LogUniformPrior("β", 0.1, 2.0),
    LogUniformPrior("β2", 0.1, 5.0),
    LogUniformPrior("ϵ", 0.1, 5.0), # constrain ϵ<5
    DependentLogUniformPrior("rs", "rc", 1400.0),
    LogUniformPrior("T0", 1.0, 50.0),
    LogUniformPrior("Tmin/T0", 0.01, 2.0),
    LogUniformPrior("rcool", 1.0, 250.0),
    LogUniformPrior("acool", 0.1, 12.0),
    LogUniformPrior("rt", 0.01, 10.0),
    UniformPrior("a", -1.0, 1.0),
    LogUniformPrior("b", 0.1, 6.0),
    LogUniformPrior("c", 0.1, 12.0),
]

# setting overdensity to 200 instead of default 500
# model(args...; kwargs...) = Model_NFW(args...; Δ=200, kwargs...)

sample(
    data,
    (1.0u"keV", 4.0u"keV"),
    Model_Vikhlinin2006,
    priors_v2006,
    0.022e22u"cm^-2",
    0.1,
    (1400, 3500),
    (1400, 3500);
    bin_size=50,
    centre_radius=0,
    use_interpolation=false,
    # use_stepsampler=true,
    # ultranest_run_args=(region_class=mlfriends.SimpleRegion,)
)

