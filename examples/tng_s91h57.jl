using BayesJ
using Unitful, UnitfulAstro, DimensionfulAngles
using Cosmology

include("metallicity.jl")

# Load data
data = FITSData(
    "../data/tng_s91_h57_3/y/chandra_acisi_cy0_obs_evts.fits",
    "../data/tng_s91_h57_3/bg/chandra_acisi_cy0_bg_evts.fits",
    "../data/tng_s91_h57_3/response_files/acisi/acisi_aimpt_cy0.arf",
    "../data/tng_s91_h57_3/response_files/acisi/acisi_aimpt_cy0.rmf",
    0.492u"arcsecondᵃ"
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

redshift = 0.1
centre_radius = 7
bin_size = 103
centre_radius_kpc = ustrip(u"kpc", centre_radius * bin_size * ustrip(u"radᵃ", data.pixel_edge_angle) * angular_diameter_dist(BayesJ.cosmo, redshift))
BayesJ.BayesJ.@mpiinfo "Fixed inner radius from centre exclusion." centre_radius_kpc

# sample prior set for piecewise fit
priors_piecewise = [
    DeltaPrior("x0", 80.534), DeltaPrior("y0", 67.006), # from centering finding fit
    DeltaPrior("r500", 885.7), DeltaPrior("ρ500", 2.195e-28), DeltaPrior("T500", 1.70),
    DeltaPrior("r200", 1338.2), UniformPrior("ρ200", 1.e-29, 1.e-28), UniformPrior("T200", 0.0, 5.0),
    DeltaPrior("r_edge", 2000.0), UniformPrior("ρ_edge", 1.e-29, 1.e-28), UniformPrior("T_edge", 0.0, 5.0)
]

sample(
    data,
    (0.1u"keV", 7.0u"keV"),
    Model_Piecewise,
    priors_piecewise,
    0.022e22u"cm^-2",
    redshift,
    (1340, 3400),
    (1340, 3400);
    bin_size=bin_size,
    centre_radius=centre_radius,
    abundances=abundances,
    use_interpolation=false,
    use_stepsampler=false,
    log_dir="../logs/s91h57_edge",
    resume="subfolder",
    ultranest_run_args=(max_num_improvement_loops=3, min_num_live_points=100, show_status=true),
)
