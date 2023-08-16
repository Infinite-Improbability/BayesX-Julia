include("../src/BayesJ.jl")
using .BayesJ
using Unitful, UnitfulAstro, DimensionfulAngles
using Test

function test_nfw()
    mass = 5e14u"Msun"
    fg = 0.13
    gnfw = [1.0510, 5.4905, 0.3081, 1.177] # Universal values from Arnaud 2010
    z = 0.5

    t, d = Model_NFW(mass, fg, gnfw..., z=z)

    test_model(t, d)

    tunitless, dunitless = Model_NFW(ustrip(u"Msun", mass), fg, gnfw..., z=z)

    @testset "Unitless model call matches unitful call" begin
        for r in (0.001:10)u"Mpc"
            @test t(r) == tunitless(r)
            @test d(r) == dunitless(r)
        end
    end

    return t, d
end

function test_model(temperature, density)
    @testset "Temperature and density functions" begin

        @testset "Testing at assorted radii" begin
            for r in (0.001:10)u"Mpc"
                t = temperature(r)
                d = density(r)

                @test dimension(t) == dimension(1u"keV")
                @test dimension(d) == dimension(1u"kg/m^3")

                @test isa(ustrip(u"keV", t), Real)
                @test isa(ustrip(u"kg/m^3", d), Real)

                @test isfinite(t)
                @test isfinite(d)

                @test t >= 0u"keV"
                @test d >= 0u"kg/m^3"

                @test t == temperature(-r)
                @test d == density(-r)

                @test t == temperature(eval(r))
                @test d == density(eval(r))

            end
        end

        @testset "Surface Brightness" begin
            z = 0.5
            emission_model = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                (0.1:7)u"keV",
                z,
                use_interpolation=false
            )

            s_distant = BayesJ.surface_brightness(
                100u"Mpc",
                temperature,
                density,
                z,
                Quantity(Inf, u"Mpc"),
                emission_model,
                0.492u"arcsecondᵃ"
            )

            @test all(i -> isapprox(i, 0u"cm^-2/s", atol=1e-10u"cm^-2/s"), s_distant)

            s_zero = BayesJ.surface_brightness(
                0u"Mpc",
                temperature,
                density,
                z,
                Quantity(Inf, u"Mpc"),
                emission_model,
                0.492u"arcsecondᵃ"
            )
            @test all(isfinite, s_zero)
        end

        # @testset "Behaviour at high radius" begin
        #     @test uconvert(u"keV", temperature(1000u"Mpc")) < 1e-10u"keV"
        #     @test uconvert(u"kg/m^3", density(1000u"Mpc")) ≈ 0u"kg/m^3"
        # end
    end

end

test_nfw();