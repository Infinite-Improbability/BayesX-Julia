using BayesJ
using BenchmarkTools
using Unitful, UnitfulAstro, DimensionfulAngles
using Profile

energy_bins = range(0u"keV", 13u"keV", length=1000)

T = 1.0u"keV"
r = 0.5u"kpc"
ρ = 1.0e-24u"g/cm^3"

z = 0.1
shape = (100, 100)

pixel_edge_angle = 0.0492u"arcsecondᵃ"
exposure_time = 3.0u"s"
response_function = rand(Float64, (100, length(energy_bins) - 1)) * 1u"cm^2"
centre_radius = 0
integration_limit = 10.0u"kpc"
projected_radius = 0.3u"kpc"

emission_model = BayesJ.prepare_model_mekal(
    2.2e20u"cm^-2",
    energy_bins,
    0.1,
)

# emission_model = BayesJ.prepare_model_mekal_interpolation(
#     2.2e20u"cm^-2",
#     energy_bins,
#     0.1,
# )


flux = zeros(Float32, length(energy_bins) - 1)

println("Emissions Model")
display(@benchmark emission_model($flux, $T, $ρ))

# CONSTANT MODEL
println("\n\nConstant Model")
temperature_cnst, density_cnst = Model_Constant(r, T, ρ)
display(
    @benchmark begin
        BayesJ.surface_brightness(
            $projected_radius,
            $temperature_cnst,
            $density_cnst,
            $z,
            $integration_limit,
            $emission_model,
            $pixel_edge_angle,
            $flux
        )
    end
)
display(
    @benchmark begin
        BayesJ.make_observation(
            $temperature_cnst,
            $density_cnst,
            $z,
            $shape,
            $pixel_edge_angle,
            $emission_model,
            $exposure_time,
            $response_function,
            $(0u"arcsecondᵃ", 0u"arcsecondᵃ"),
            $centre_radius,
            limit=$integration_limit,
        )
    end
)

# PIECEWISE MODEL
println("\n\nPiecewise Model")
temperature_pw, density_pw = Model_Piecewise(
    0.0, 1.e-23, 5.0,
    10.0, 1.e-24, 4.0,
    100.0, 1.e-25, 3.0,
    1000.0, 0.0, 0.0
)
display(
    @benchmark begin
        BayesJ.surface_brightness(
            $projected_radius,
            $temperature_pw,
            $density_pw,
            $z,
            $integration_limit,
            $emission_model,
            $pixel_edge_angle,
            $flux
        )
    end
)
display(
    @benchmark begin
        BayesJ.make_observation(
            $temperature_pw,
            $density_pw,
            $z,
            $shape,
            $pixel_edge_angle,
            $emission_model,
            $exposure_time,
            $response_function,
            $(0u"arcsecondᵃ", 0u"arcsecondᵃ"),
            $centre_radius,
            limit=$integration_limit,
        )
    end
)