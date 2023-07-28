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

energy_range = (0.3:0.01:0.46)u"keV"
exposure_time = 300e3u"s"
pixel_size = 0.492u"arcsecondᵃ"
redshift = 0.5

mass = 7e14u"Msun"
fg = 0.13

response = load_response(data, energy_range)
em = prepare_model_mekal(2.2e20u"cm^-2", energy_range)
model = Model_NFW_GNFW(mass, fg, 1.062, 5.4807, 0.3292, 1.156, redshift, [64, 64], pixel_size, em, exposure_time, response)

bg_rate = 8.4e-6u"cm^-2/arcminuteᵃ^2/s"
avg_eff_area = 250u"cm^2"
n_channels = size(model, 1)
bg_count = bg_rate * exposure_time * avg_eff_area * pixel_size^2 / n_channels
bg = fill(bg_count, size(model))

s = sample(
    round.(Int, model + bg),
    round.(Int, ustrip.(Float64, NoUnits, bg)),
    response,
    make_cube_transform(LogUniformPrior(1e13, 1e15), UniformPrior(0.01, 1.0)),
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
