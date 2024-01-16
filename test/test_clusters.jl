# include("../src/BayesJ.jl")
# using .BayesJ
using Unitful, UnitfulAstro, DimensionfulAngles

logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))
logrange(x1::Unitful.Length, x2::Unitful.Length, n) = (1u"kpc" * 10^y for y in range(log10(ustrip(u"kpc", x1)), log10(ustrip(u"kpc", x2)), length=n))

function test_model(temperature, density)
    @testset "Temperature and density functions" begin

        @testset "Testing at assorted radii" begin
            for r in logrange(1u"pc", 1u"Mpc", 100)
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
            energy_range = (0.1:7)u"keV"
            emission_model = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_range,
                z,
            )
            flux = zeros(Float32, length(energy_range) - 1)

            @testset "Zero at distance" begin
                s_distant = BayesJ.surface_brightness(
                    100u"Mpc",
                    temperature,
                    density,
                    z,
                    Quantity(1000.0, u"Mpc"),
                    emission_model,
                    0.492u"arcsecondᵃ",
                    flux
                )


                @test all(i -> isapprox(i, 0u"cm^-2/s", atol=1e-10u"cm^-2/s"), s_distant)
                if !(all(i -> isapprox(i, 0u"cm^-2/s", atol=1e-10u"cm^-2/s"), s_distant))
                    display(s_distant)
                end
            end

            @testset "Finite at centre" begin
                s_zero = BayesJ.surface_brightness(
                    0u"Mpc",
                    temperature,
                    density,
                    z,
                    Quantity(10.0, u"Mpc"),
                    emission_model,
                    0.492u"arcsecondᵃ",
                    flux
                )

                @test all(isfinite, s_zero)
            end

        end

        # @testset "Behaviour at high radius" begin
        #     @test uconvert(u"keV", temperature(1000u"Mpc")) < 1e-10u"keV"
        #     @test uconvert(u"kg/m^3", density(1000u"Mpc")) ≈ 0u"kg/m^3"
        # end
    end

end

function test_nfw()
    @testset "NFW" begin
        mass = 5e14u"Msun"
        fg = 0.13
        c_dm = 3.0
        gnfw = [1.0510, 5.4905, 0.3081, 1.177] # Universal values from Arnaud 2010
        z = 0.5

        t, d = Model_NFW(mass, fg, c_dm, gnfw..., z=z)

        test_model(t, d)

        @testset "Unitless model call matches unitful call" begin
            tunitless, dunitless = Model_NFW(ustrip(u"Msun", mass), fg, c_dm, gnfw..., z=z)
            for r in logrange(1u"pc", 1u"Mpc", 100)
                @test t(r) == tunitless(r)
                @test d(r) == dunitless(r)
            end
        end
    end
end

function test_einasto()
    @testset "Einasto" begin
        mass = 5e14u"Msun"
        fg = 0.13
        c_dm = 3.0
        gnfw = [1.0510, 5.4905, 0.3081, 1.177] # Universal values from Arnaud 2010
        z = 0.5

        for α in 0.6:0.1:1.9
            @testset "α=$α" begin
                t, d = Model_Einasto(mass, fg, c_dm, α, gnfw..., z=z)

                test_model(t, d)

                @testset "Unitless model call matches unitful call" begin
                    tunitless, dunitless = Model_Einasto(ustrip(u"Msun", mass), fg, c_dm, α, gnfw..., z=z)
                    for r in logrange(1u"pc", 1u"Mpc", 100)
                        @test t(r) == tunitless(r)
                        @test d(r) == dunitless(r)
                    end
                end
            end
        end

        @test_throws BayesJ.PriorError Model_Einasto(mass, fg, c_dm, 2.1, gnfw..., z=z)
        @test_throws BayesJ.PriorError Model_Einasto(mass, fg, c_dm, 0.01, gnfw..., z=z)
    end
end

function test_vikhlinin2006()
    @testset "Vikhlinin2006" begin
        function strip_units(i)
            if dimension(i) == Unitful.𝐋^-3
                return ustrip(u"cm^-3", i)
            elseif isa(i, Unitful.Length)
                return ustrip(u"kpc", i)
            elseif isa(i, Unitful.Energy)
                return ustrip(u"keV", i)
            else
                return i
            end
        end

        params = (
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
            10.0,
            1.e-9,
        ) # Universal values from Arnaud 2010
        unitless_params = Tuple(strip_units(i) for i in params)
        z = 0.5

        t, d = Model_Vikhlinin2006(params..., z=z)
        # display(t(100u"kpc"))
        # display(d(100u"kpc"))

        test_model(t, d)

        @testset "Unitless model call matches unitful call" begin
            tunitless, dunitless = Model_Vikhlinin2006(unitless_params..., z=z)
            for r in logrange(1u"pc", 1u"Mpc", 100)
                @test t(r) == tunitless(r)
                @test d(r) == dunitless(r)
            end
        end
    end

end

@testset "Cluster Models" begin
    test_nfw()
    test_einasto()
    test_vikhlinin2006()
end