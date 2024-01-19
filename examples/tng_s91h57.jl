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

priors_piecewise = [
    DeltaPrior("x0", 80.534), DeltaPrior("y0", 67.006), # from centering finding fit
    DeltaPrior("r0", 0.0), UniformPrior("ρ0", 0.0, 1e-22), UniformPrior("T0", 0.0, 10.0),
    UniformPrior("r1", 0.0, 1000.0), UniformPrior("ρ1", 0.0, 1e-22), UniformPrior("T1", 0.0, 10.0),
    DependentUniformPrior("r2", "r1", 1000.0), UniformPrior("ρ2", 0.0, 1e-22), UniformPrior("T2", 0.0, 10.0),
    DependentUniformPrior("r3", "r2", 1000.0), UniformPrior("ρ3", 0.0, 1e-22), UniformPrior("T3", 0.0, 10.0),
    DependentUniformPrior("r4", "r3", 1000.0), DeltaPrior("ρ4", 0.0), DeltaPrior("T4", 0.0),
]

sample(
    data,
    (0.7u"keV", 7.0u"keV"),
    Model_Piecewise,
    priors_piecewise,
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
