using Plots
using Profile
using BenchmarkTools
using DotEnv
using Unitful, DimensionfulAngles
using PoissonRandom
# using FortranFiles

# We use include instead of importing the module for easy access to internal methods
include("../src/run.jl")

# f = FortranFile("../BayesX/log.dat", "r")
# fm = reshape(read(f, (Float64, 64 * 64 * 32)), (32, 64, 64))

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
pixel_size = 0.492e1u"arcsecondᵃ"
redshift = 0.164
shape = [32, 32]
gnfw = [1.0510, 5.4905, 0.3081, 1.177] # Using universal values from Arnaud 2010
mass = 5e14u"Msun"
fg = 0.13

priors = Model_NFW_GNFW_Priors(
    UniformPrior(1e14, 1e15),
    UniformPrior(0.08, 0.2),
    DeltaPrior(gnfw[1]),
    DeltaPrior(gnfw[2]),
    DeltaPrior(gnfw[3]),
    DeltaPrior(gnfw[4]),
    DeltaPrior(0),
    DeltaPrior(0)
)

response = load_response(data, energy_range)
em = prepare_model_mekal(
    2.2e20u"cm^-2",
    energy_range,
    redshift,
    use_interpolation=true,
    temperatures=(1e-30:0.05:9.0)u"keV",
    hydrogen_densities=(1e-30:0.005:1.0)u"cm^-3"
)
model = Model_NFW_GNFW(
    MT_200=mass,
    fg_200=fg,
    α=gnfw[1],
    β=gnfw[2],
    γ=gnfw[3],
    c_500_GNFW=gnfw[4],
    z=redshift,
    shape=shape,
    pixel_edge_angle=pixel_size,
    emission_model=em,
    exposure_time=exposure_time,
    response_function=response,
    centre_coordinates=(0u"arcsecondᵃ", 0u"arcsecondᵃ"),
    centre_radius=0
)
replace!(model, NaN => 0.0)

@assert all(isfinite, model)

noisy = pois_rand.(model)

# @profview Model_NFW_GNFW(
#     mass,
#     fg,
#     gnfw...,
#     redshift,
#     shape,
#     pixel_size,
#     em,
#     exposure_time,128
#     (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
#     0
# )
# @profview_allocs Model_NFW_GNFW(
#     mass,
#     fg,
#     gnfw...,
#     redshift,
#     shape,
#     pixel_size,
#     em,
#     exposure_time,
#     response,
#     (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
#     0
# )
# display(
#     @benchmark Model_NFW_GNFW(
#         mass,
#         fg,
#         gnfw...,
#         redshift,
#         shape,
#         pixel_size,
#         em,
#         exposure_time,
#         response,
#         (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
#         0
#     )
# )

bg_rate = 8.4e-6u"cm^-2/arcminuteᵃ^2/s";
avg_eff_area = 250u"cm^2";
n_channels = size(model, 1);
bg_count = bg_rate * exposure_time * avg_eff_area *
           pixel_size^2 / n_channels;
bg = fill(upreferred(bg_count), size(model));

# display(
#     heatmap(
#         dropdims(sum(model + bg, dims=1), dims=1),
#         title="Model"
#     )
# )
# display(
#     heatmap(
#         dropdims(sum(noisy + bg, dims=1), dims=1),
#         title="Noisy"
#     )
# )

s = sample(
    round.(Int, noisy + bg),
    round.(Int, bg),
    response,
    generate_transform(priors),
    exposure_time,
    exposure_time,
    redshift,
    emission_model=em,
    pixel_edge_angle=pixel_size,
    background_rate=bg_rate,
    average_effective_area=avg_eff_area,
    centre_radius=12
)

posterior = s[2]["posterior"]
@assert all(posterior["errlo"][1:2] .< [ustrip(mass), fg] .< posterior["errup"][1:2]) display(posterior)