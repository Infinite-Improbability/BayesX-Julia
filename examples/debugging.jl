using Plots
using Profile
using BenchmarkTools
using DotEnv
using Unitful, DimensionfulAngles
using PoissonRandom
using FortranFiles

# We use include instead of importing the module for easy access to internal methods
include("../src/run.jl")

f = FortranFile("../BayesX/log.dat", "r")
fm = reshape(read(f, (Float64, 64 * 64 * 32)), (32, 64, 64))

DotEnv.config()

data = PlaintextData(
    "data/simtestdata/data64by64by32.txt",
    "data/simtestdata/BG64by64by32.txt",
    "data/simtestdata/ARF_32by1.txt",
    "data/simtestdata/RMF_32by32.txt",
    (32, 64, 64),
    (32, 32),
    300e3u"s",
    300e3u"s",
    0.492u"arcsecondᵃ",
    1u"cm^2"
)

energy_range = range(0.3u"keV", 7.0u"keV", 33)
exposure_time = 300e3u"s"
pixel_size = 0.492u"arcsecondᵃ"
redshift = 0.164

mass = 3e14u"Msun"
fg = 0.13

response = load_response(data, energy_range)
em = prepare_model_mekal(2.2e20u"cm^-2", energy_range, temperatures=(0.001:0.005:20)u"keV", hydrogen_densities=(1e-10:0.005:10)u"cm^-3", use_interpolation=true)
model = Model_NFW_GNFW(mass, fg, 1.062, 5.4807, 0.3292, 1.156, redshift, [64, 64], pixel_size, em, exposure_time, response)
noisy = pois_rand.(model)

# @profview Model_NFW_GNFW(mass, fg, 1.062, 5.4807, 0.3292, 1.156, redshift, [24, 24], pixel_size, em, exposure_time, response)
# @profview_allocs Model_NFW_GNFW(mass, fg, 1.062, 5.4807, 0.3292, 1.156, redshift, [24, 24], pixel_size, em, exposure_time, response)
display(@benchmark Model_NFW_GNFW(mass, fg, 1.062, 5.4807, 0.3292, 1.156, redshift, [24, 24], pixel_size, em, exposure_time, response))

bg_rate = 8.4e-6u"cm^-2/arcminuteᵃ^2/s";
avg_eff_area = 250u"cm^2";
n_channels = size(model, 1);
bg_count = bg_rate * exposure_time * avg_eff_area * pixel_size^2 / n_channels;
bg = fill(upreferred(bg_count), size(model));

# display(heatmap(dropdims(sum(model + bg, dims=1), dims=1), title="Model"))
# display(heatmap(dropdims(sum(noisy + bg, dims=1), dims=1), title="Noisy"))

# s = sample(
#     round.(Int, noisy + bg),
#     round.(Int, bg),
#     response,
#     make_cube_transform(UniformPrior(1e14, 1e15), UniformPrior(0.08, 0.2)),
#     exposure_time,
#     exposure_time,
#     redshift,
#     emission_model=em,
#     pixel_edge_angle=pixel_size,
#     background_rate=bg_rate,
#     average_effective_area=avg_eff_area
# )

# posterior = s[2]["posterior"]
# @assert all(posterior["errlo"] .< [ustrip(mass), fg] .< posterior["errup"]) display(posterior)

# The ongoing problem seems to be that count rates are just too low to get
# information with realistic exposure times. But clearly BayesX lacks this
# problem.
# Adding noise may have fixed this