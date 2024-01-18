using BayesJ
using Unitful, DimensionfulAngles

include("metallicity.jl")

data = FITSData(
    "data/tng/tng_s91_h57_y_obs_evt.fits",
    "data/tng/tng_s91_h57_y_bg.fits",
    "data/tng/acisi_aimpt_cy0.arf",
    "data/tng/acisi_aimpt_cy0.rmf",
    0.492u"arcsecondᵃ"
)

gas_metal_fractions = [
    0.7543419003486633,
    0.24289754033088684,
    0.0002910351031459868,
    9.24545165617019e-05,
    0.001410270226188004,
    0.0004578478983603418,
    0.00013437608140520751,
    0.00012510872329585254,
    0.00016363433678634465,
    8.585586328990757e-05
]
gas_metals = [
    "H", "He", "C", "N", "O", "Ne", "Mg", "Si", "Fe", "Other"
]
abundances = convert_to_anders(gas_metals, gas_metal_fractions)

priors_v2006 = [
    DeltaPrior("x0", 80.534), # from centre finding fit
    DeltaPrior("y0", 67.006), # from centering finding fit
    UniformPrior("n0", 0.1e-3, 40.0e-3),
    UniformPrior("n02", 0.001e-2, 6.0e-1),
    DependentLogUniformPrior("rc", "rc2", 600.0),
    UniformPrior("rc2", 1.0, 100.0),
    UniformPrior("α", 0.1, 3.0),
    UniformPrior("β", 0.1, 2.0),
    UniformPrior("β2", 0.1, 5.0),
    UniformPrior("ϵ", 0.1, 5.0), # constrain ϵ<5
    DependentLogUniformPrior("rs", "rc", 1400.0),
    UniformPrior("T0", 1.0, 20.0),
    UniformPrior("Tmin/T0", 0.01, 1.0),
    UniformPrior("rcool", 1.0, 250.0),
    UniformPrior("acool", 0.1, 10.0),
    UniformPrior("rt", 0.01, 10.0),
    UniformPrior("a", -1.0, 1.0),
    UniformPrior("b", 0.1, 6.0),
    UniformPrior("c", 0.1, 10.0),
    LogUniformPrior("d", 1.0e-12, 1.0),
]

sample(
    data,
    (0.7u"keV", 7.0u"keV"),
    Model_Vikhlinin2006,
    priors_v2006,
    0.022e22u"cm^-2",
    0.1,
    (1340, 3400),
    (1340, 3400);
    bin_size=40,
    centre_radius=0,
    abundances=abundances,
    use_interpolation=false,
    use_stepsampler=false,
    log_dir="logs/s91h57_vikh",
    ultranest_run_args=(max_num_improvement_loops=3, min_num_live_points=400, show_status=false),
)
