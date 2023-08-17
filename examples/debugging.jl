using Plots
using Profile
using BenchmarkTools
using DotEnv
using Unitful, DimensionfulAngles
using PoissonRandom

ENV["JULIA_DEBUG"] = "Main"

# We use include instead of importing the module for easy access to internal methods
include("../src/run.jl")
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

response = load_response(data, energy_range)
em = prepare_model_mekal(
    2.2e20u"cm^-2",
    energy_range,
    redshift,
    use_interpolation=true,
    temperatures=(1e-30:0.05:9.0)u"keV",
    hydrogen_densities=(1e-30:0.005:3.0)u"cm^-3"
)
# cluster = Model_NFW(
#     mass,
#     fg,
#     gnfw...,
#     redshift,ulltext/55351.text.html
# )
cluster = Model_Vikhlinin2006(
    4.705e-3u"cm^-3",
    0.247e-1u"cm^-3",
    94.6u"kpc",
    75.83u"kpc",
    0.916,
    0.526,
    3.607,
    4.943,
    1239.9u"kpc",
    3.61u"keV",
    0.27,
    57u"kpc",
    3.88,
    1.42u"Mpc",
    0.12,
    5.00,
    10.0
)

@info "Profiles generated"
obs = make_observation(
    cluster...,
    redshift,
    shape,
    pixel_size,
    em,
    exposure_time,
    response,
    (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
    1
)

replace!(obs, NaN => 0.0)
@assert all(isfinite, obs)
@info "Observation generated"
ENV["JULIA_DEBUG"] = ""

@info "Running tests"

function test_func()
    cluster = Model_Einasto(
        mass,
        fg,
        0.5,
        gnfw...,
        z=0.5
    )
    make_observation(
        cluster...,
        redshift,
        shape,
        pixel_size,
        em,
        exposure_time,
        response,
        (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
        1
    )
end

display(
    @benchmark test_func()
)
@profview test_func()

@info "Preparing data for sampling"

noisy = pois_rand.(obs)
bg_rate = 8.4e-6u"cm^-2/arcminuteᵃ^2/s";
avg_eff_area = 250u"cm^2";
n_channels = size(obs, 1);
bg_count = bg_rate * exposure_time * avg_eff_area *
           pixel_size^2 / n_channels;
bg = fill(upreferred(bg_count), size(obs));

display(heatmap(dropdims(sum(obs + bg, dims=1), dims=1), title="Model"))
display(heatmap(dropdims(sum(noisy + bg, dims=1), dims=1), title="Noisy"))

T(r) = uconvert(u"keV", cluster[1](r))
n(r) = uconvert(u"cm^-3", cluster[2](r))

# s = sample(
#     round.(Int, noisy + bg),
#     round.(Int, bg),
#     response,
#     make_cube_transform(
#         UniformPrior("x0", -1.0, 1.0),
#         UniformPrior("y0", -1.0, 1.0),
#         UniformPrior("MT_200", 1e14, 1e15),
#         UniformPrior("fg_200", 0.08, 0.2),
#     ),
#     exposure_time,
#     exposure_time,
#     redshift,
#     prior_names=["x0", "y0", "MT_200", "fg_200"],
#     cluster_model=Model_NFW,
#     emission_model=em,
#     pixel_edge_angle=pixel_size,
#     background_rate=bg_rate,
#     average_effective_area=avg_eff_area,
#     centre_radius=1
# )

# posterior = s[2]["posterior"]
# @assert all(posterior["errlo"][1:2] .< [ustrip(mass), fg] .< posterior["errup"][1:2]) display(posterior)