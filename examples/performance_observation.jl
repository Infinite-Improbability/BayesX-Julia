using Profile, BenchmarkTools
using Unitful, DimensionfulAngles, UnitfulAstro

using BayesJ

@info "Start"

data = PlaintextData(
    "../data/simtestdata/data64by64by32.txt",
    "../data/simtestdata/BG64by64by32.txt",
    "../data/simtestdata/ARF_32by1.txt",
    "../data/simtestdata/RMF_32by32.txt",
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
cΔ = 3.0

response = BayesJ.load_response(data, energy_range)
em = BayesJ.prepare_model_mekal(
    2.2e20u"cm^-2",
    energy_range,
    redshift,
    use_interpolation=false,
    temperatures=(1e-30:0.05:9.0)u"keV",
    hydrogen_densities=(1e-30:0.005:3.0)u"cm^-3",
)

cluster_nfw = Model_NFW(
    mass,
    fg,
    cΔ,
    gnfw...,
    z=0.5
)

cluster_einasto = Model_Einasto(
    mass,
    fg,
    cΔ,
    0.5,
    gnfw...,
    z=0.5
)

cluster_vk2006 = Model_Vikhlinin2006(
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

limFinite = Quantity(5.0, u"Mpc")
limInfinite = Quantity(Inf, u"Mpc")

@info "Infinite"
display(
    @benchmark BayesJ.make_observation(
        $cluster_vk2006...,
        $redshift,
        $shape,
        $pixel_size,
        $em,
        $exposure_time,
        $response,
        (0, 0),
        1,
        limit=$limInfinite
    )
)

@info "Finite"
display(
    @benchmark BayesJ.make_observation(
        $cluster_vk2006...,
        $redshift,
        $shape,
        $pixel_size,
        $em,
        $exposure_time,
        $response,
        (0, 0),
        1,
        limit=$limFinite
    )
)

# @profview BayesJ.make_observation(
#     cluster...,
#     redshift,
#     shape,
#     pixel_size,
#     em,
#     exposure_time,
#     response,
#     (0, 0),
#     1,
#     limit=limInf
# )
