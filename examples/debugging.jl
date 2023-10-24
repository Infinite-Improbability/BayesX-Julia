using Plots
using Profile
using BenchmarkTools
using Unitful, DimensionfulAngles, UnitfulAstro
using PoissonRandom

# include("../src/BayesJ.jl")

using BayesJ

ENV["JULIA_DEBUG"] = "BayesJ"

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
shape = [256, 256]
gnfw = [1.0510, 5.4905, 0.3081, 1.177] # Using universal values from Arnaud 2010
mass = 5e14u"Msun"
fg = 0.13

response = BayesJ.load_response(data, energy_range)
em = BayesJ.prepare_model_mekal(
    2.2e20u"cm^-2",
    energy_range,
    redshift,
    use_interpolation=false,
    temperatures=(1e-30:0.05:9.0)u"keV",
    hydrogen_densities=(1e-30:0.005:3.0)u"cm^-3",
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
obs = BayesJ.make_observation(
    cluster...,
    redshift,
    shape,
    pixel_size,
    em,
    exposure_time,
    response,
    (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
    0
)

replace!(obs, NaN => 0.0)  # does this need handling for missing?
@assert all(isfinite, obs)
@info "Observation generated"
ENV["JULIA_DEBUG"] = ""

@info "Running tests"

function test_func()
    cluster = Model_Einasto(
        mass,
        fg,
        3.0,
        0.5,
        gnfw...,
        z=0.5
    )
    BayesJ.make_observation(
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
