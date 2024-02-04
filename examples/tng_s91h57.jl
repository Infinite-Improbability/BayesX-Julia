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
    0.322u"arcsecondᵃ"
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
bin_size = 128
centre_radius_kpc = ustrip(u"kpc", centre_radius * bin_size * ustrip(u"radᵃ", data.pixel_edge_angle) * angular_diameter_dist(BayesJ.cosmo, redshift))
inner_radius = centre_radius_kpc * 1.1
r1 = (inner_radius + 1200.0) / 2
BayesJ.BayesJ.@mpiinfo "Fixed inner radius from centre exclusion." centre_radius_kpc inner_radius r1

# sample prior set for piecewise fit
priors_piecewise = [
    DeltaPrior("x0", 0.0), DeltaPrior("y0", 0.0),
    DeltaPrior("r0", inner_radius), DeltaPrior("ρ0", 5.24e-28), DeltaPrior("T0", 2.074),
    DeltaPrior("r1", r1), UniformPrior("ρ1", 1.e-28, 1.e-27), UniformPrior("T1", 1.407, 2.074),
    DeltaPrior("r2", 1200.0), DeltaPrior("ρ2", 1.137e-28), DeltaPrior("T2", 1.407),
]

sample(
    data,
    (0.05u"keV", 2.5u"keV"),
    Model_Piecewise,
    priors_piecewise,
    0.022e22u"cm^-2",
    redshift,
    (2048, 6144),
    (2048, 6144);
    bin_size=bin_size,
    centre_radius=centre_radius,
    abundances=abundances,
    use_interpolation=false,
    use_stepsampler=false,
    log_dir="../logs/s91h57_edge",
    resume="subfolder",
    ultranest_run_args=(max_num_improvement_loops=3, min_num_live_points=100, show_status=true),
)
