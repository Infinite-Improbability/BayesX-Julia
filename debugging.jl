using Plots
using Profile
using BenchmarkTools
using DotEnv
using Unitful, DimensionfulAngles

include("src/run.jl")

# f = FortranFile("../BayesX/log.dat")
# fm = reshape(read(f, (Float64, 64 * 64 * 32)), (32, 64, 64))

DotEnv.config()

data = FITSData(
    "$(ENV["DEBUG_DATA_PATH"])/acisf04361_repro_evt2.fits",
    "$(ENV["DEBUG_DATA_PATH"])/4361_blanksky.fits",
    "$(ENV["DEBUG_DATA_PATH"])/specx/specx.arf",
    "$(ENV["DEBUG_DATA_PATH"])/specx/specx.rmf",
    0.492u"arcsecondᵃ"
)

energy_range = (0.3:0.01:3.0)u"keV"
exposure_time = 300e6u"s" # inflated to crazy lengths
pixel_size = 0.492u"arcsecondᵃ"
redshift = 0.164

mass = 3e14u"Msun"
fg = 0.13

response = load_response(data, energy_range)
em = prepare_model_mekal(2.2e20u"cm^-2", energy_range)
model = Model_NFW_GNFW(mass, fg, 1.062, 5.4807, 0.3292, 1.156, redshift, [32, 32], pixel_size, em, exposure_time, response)

# @profview Model_NFW_GNFW(mass, fg, 1.062, 5.4807, 0.3292, 1.156, redshift, [24, 24], pixel_size, em, exposure_time, response)
# @profview_allocs Model_NFW_GNFW(mass, fg, 1.062, 5.4807, 0.3292, 1.156, redshift, [24, 24], pixel_size, em, exposure_time, response)
# @benchmark Model_NFW_GNFW(mass, fg, 1.062, 5.4807, 0.3292, 1.156, redshift, [24, 24], pixel_size, em, exposure_time, response)

bg_rate = 8.4e-6u"cm^-2/arcminuteᵃ^2/s"
avg_eff_area = 250u"cm^2"
n_channels = size(model, 1)
bg_count = bg_rate * exposure_time * avg_eff_area * pixel_size^2 / n_channels
bg = fill(upreferred(bg_count), size(model))

# TODO: The noise! I forgot the noise

s = sample(
    round.(Int, model + bg),
    round.(Int, bg),
    response,
    make_cube_transform(UniformPrior(1e14, 1e15), UniformPrior(0.08, 0.2)),
    exposure_time,
    exposure_time,
    redshift,
    emission_model=em,
    pixel_edge_angle=pixel_size,
    background_rate=bg_rate,
    average_effective_area=avg_eff_area
)

posterior = s[2]["posterior"]
@assert all(posterior["errlo"] .< [ustrip(mass), fg] .< posterior["errup"]) display(posterior)

# The ongoing problem seems to be that count rates are just too low to get
# information with realistic exposure times. But clearly BayesX lacks this
# problem.
# Also need to add noise.