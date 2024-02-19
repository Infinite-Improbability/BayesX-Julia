using BayesJ
using Unitful, UnitfulAstro, DimensionfulAngles
using Cosmology

include("metallicity.jl")

# Load data
const data = FITSData(
    "../data/tng_s91_h57_4/y/chandra_acisi_cy0_obs_evts.fits",
    "../data/tng_s91_h57_4/bg/chandra_acisi_cy0_bg_evts.fits",
    "../data/tng_s91_h57_4/response_files/acisi/acisi_aimpt_cy0.arf",
    "../data/tng_s91_h57_4/response_files/acisi/acisi_aimpt_cy0.rmf",
    0.492u"arcsecondᵃ"
)

# Load abundances from tng metadata
const gas_metal_fractions = [
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
const gas_metals = [
    "H", "He", "C", "N", "O", "Ne", "Mg", "Si", "Fe", "Other"
]
const abundances = convert_to_anders(gas_metals, gas_metal_fractions)

const c200 = 6.077
const r200 = 1338.216892862
const rs = r200 / c200
const r500 = 885.727387822
const c500 = r500 / rs

const priors_nfw = [
    NormalPrior("x0", 0.0, 25.0), NormalPrior("y0", 0.0, 25.0),
    DeltaPrior("MT_500", 20445.219658806e10),
    DeltaPrior("fg_500", 0.13),
    DeltaPrior("c_500", c500),
    DeltaPrior("a", 1.0510),
    DeltaPrior("b", 5.4905),
    DeltaPrior("c", 0.3081),
    DeltaPrior("c_500_GNFW", 1.177)
]

sample(
    data,
    (0.7u"keV", 7.0u"keV"),
    Model_NFW,
    priors_nfw,
    0.022e22u"cm^-2",
    0.1,
    (1340, 3340),
    (1370, 3400);
    bin_size=32,
    centre_radius=0,
    abundances=abundances,
    use_stepsampler=false,
    log_dir="../logs/s91h57/nfw",
    resume="subfolder",
    ultranest_run_args=(max_num_improvement_loops=3, min_num_live_points=100),
)
