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
    DeltaPrior("r0", 0.958), DeltaPrior("ρ0", 1.28e-26), DeltaPrior("T0", 5.04),
    DeltaPrior("r1", 250.0), UniformPrior("ρ1", 1.e-27, 1.e-26), UniformPrior("T1", 2.0, 3.0),
    DeltaPrior("r2", 500.0), UniformPrior("ρ2", 1.e-28, 1e-27), UniformPrior("T2", 2.0, 3.0),
    DeltaPrior("r4", 2257), DeltaPrior("ρ4", 6.35e-29), DeltaPrior("T4", 0.27),
]

priors_nfw = [
    # UniformPrior("x0", -100.0, 100.0),
    # UniformPrior("y0", -100.0, 100.0),
    DeltaPrior("x0", 47.5),
    DeltaPrior("y0", 47.5),
    DeltaPrior("MT_200", 2.8e14),
    NormalPrior("fg_500", 0.13, 0.01),
    DeltaPrior("c_200", 6.07),
    NormalPrior("a", 1.0620, 0.06),
    NormalPrior("b", 5.4807, 1.0),
    NormalPrior("c", 0.3292, 0.02),
    NormalPrior("c_500_GNFW", 1.156, 0.02)
]

model(args...; kwargs...) = Model_NFW(args...; Δ=200, kwargs...)

sample(
    data,
    (0.7u"keV", 7.0u"keV"),
    model,
    priors_nfw,
    0.022e22u"cm^-2",
    0.1,
    (1340, 3400),
    (1340, 3400);
    bin_size=120,
    centre_radius=0,
    abundances=abundances,
    use_interpolation=false,
    use_stepsampler=false,
    log_dir="logs/s91h57_nfw/run7",
    resume="resume",
    ultranest_run_args=(max_num_improvement_loops=3, min_num_live_points=100, show_status=true),
)
