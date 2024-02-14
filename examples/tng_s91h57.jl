using BayesJ
using Unitful, UnitfulAstro, DimensionfulAngles
using Cosmology

include("metallicity.jl")

# Load data
data = FITSData(
    "../data/tng_s91_h57_4/y/lynx_hdxi_obs_evts.fits",
    "../data/tng_s91_h57_4/bg/lynx_hdxi_bg_evts.fits",
    "../data/tng_s91_h57_4/response_files/lynx_hdxi/xrs_hdxi_3x10.arf",
    "../data/tng_s91_h57_4/response_files/lynx_hdxi/xrs_hdxi.rmf",
    0.33u"arcsecondᵃ"
)

# Load abundances from tng metadata
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

priors_nfw = [
    NormalPrior("x0", 0.0, 25.0), NormalPrior("y0", 0.0, 25.0),
    UniformPrior("MT_200", 1.0e14, 4.0e14),
    NormalPrior("fg_200", 0.13, 0.01),
    DeltaPrior("c_200", 6.077),
    NormalPrior("a", 1.0510, 0.06),
    NormalPrior("b", 5.4905, 1.0),
    NormalPrior("c", 0.3081, 0.02),
    NormalPrior("c_500_GNFW", 1.177, 0.02)
]

priors_einasto = [
    NormalPrior("x0", 0.0, 25.0), NormalPrior("y0", 0.0, 25.0),
    UniformPrior("MT_200", 1.0e14, 4.0e14),
    NormalPrior("fg_200", 0.13, 0.01),
    DeltaPrior("c_200", 6.077),
    UniformPrior("n", 0.1, 10.0),
    NormalPrior("a", 1.0510, 0.06),
    NormalPrior("b", 5.4905, 1.0),
    NormalPrior("c", 0.3081, 0.02),
    NormalPrior("c_500_GNFW", 1.177, 0.02)
]

model_nfw(args...; kwargs...) = Model_NFW(args...; Δ=200, kwargs...)
model_einasto(args...; kwargs...) = Model_Einasto(args...; Δ=200, kwargs...)

sample(
    data,
    (0.7u"keV", 7.0u"keV"),
    model_nfw,
    priors_nfw,
    0.022e22u"cm^-2",
    0.1,
    (2048, 6144),
    (2048, 6144);
    bin_size=128,
    centre_radius=0,
    abundances=abundances,
    use_stepsampler=false,
    log_dir="../logs/s91h57/nfw",
    resume="subfolder",
    ultranest_run_args=(max_num_improvement_loops=3, min_num_live_points=100),
)
