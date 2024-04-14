# include("../src/BayesJ.jl")
# using .BayesJ
using Unitful, UnitfulAstro, DimensionfulAngles

logrange(x1, x2, n) = (10^y for y in range(log10(x1), log10(x2), length=n))
logrange(x1::Unitful.Length, x2::Unitful.Length, n) = (1u"kpc" * 10^y for y in range(log10(ustrip(u"kpc", x1)), log10(ustrip(u"kpc", x2)), length=n))

function test_model(cluster_model::ClusterModel)
    @testset "Temperature at various radii" begin
        for r in logrange(1u"pc", 1u"Mpc", 100)
            t = temperature(cluster_model, r)

            @test dimension(t) == dimension(1u"keV")
            @test isa(ustrip(u"keV", t), Real)
            @test isfinite(t)
            @test t >= 0u"keV"
            # @test t == temperature(-r)
            # Testing result stays consistent
            # Eval to prevent the compiler playing any tricks
            @test t == temperature(cluster_model, eval(r))

        end
    end

    @testset "Density at various radii" begin
        for r in logrange(1u"pc", 1u"Mpc", 100)
            d = density(cluster_model, r)

            @test dimension(d) == dimension(1u"kg/m^3")
            @test isa(ustrip(u"kg/m^3", d), Real)
            @test isfinite(d)
            @test d >= 0u"kg/m^3"
            # @test d == density(-r)
            @test d == density(cluster_model, eval(r))
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
                cluster_model,
                z,
                Quantity(1000.0, u"Mpc"),
                emission_model,
                0.492u"arcsecondᵃ",
                flux
            )


            @test all(i -> isapprox(i, 0u"cm^-2/s", atol=1e-10u"cm^-2/s"), s_distant)
            # if !(all(i -> isapprox(i, 0u"cm^-2/s", atol=1e-10u"cm^-2/s"), s_distant))
            #     display(s_distant)
            # end
        end

        @testset "Finite at centre" begin
            s_zero = BayesJ.surface_brightness(
                0u"Mpc",
                cluster_model,
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

function unitless_test(model::C, unitless::C) where {C<:ClusterModel}
    @testset "Unitless model call matches unitful call" begin
        for r in logrange(1u"pc", 1u"Mpc", 100)
            @test temperature(model, r) == temperature(unitless, r)
            @test density(model, r) == density(unitless, r)
        end
    end
end

function test_nfw()
    @testset "NFW" begin
        mass = 5e14u"Msun"
        fg = 0.13
        c_dm = 3.0
        gnfw = [1.0510, 5.4905, 0.3081, 1.177] # Universal values from Arnaud 2010
        z = 0.5

        model = NFWModel(mass, fg, c_dm, gnfw..., z=z)
        test_model(model)

        unitless = NFWModel(ustrip(u"Msun", mass), fg, c_dm, gnfw..., z=z)
        unitless_test(model, unitless)
    end
end

function test_einasto()
    @testset "Einasto" begin
        mass = 5e14u"Msun"
        fg = 0.13
        c_dm = 3.0
        gnfw = [1.0510, 5.4905, 0.3081, 1.177] # Universal values from Arnaud 2010
        z = 0.5

        for n in 0.3:0.5:10
            @testset "n=$n" begin
                model = EinastoModel(mass, fg, c_dm, n, gnfw..., z=z)
                test_model(model)

                unitless = EinastoModel(ustrip(u"Msun", mass), fg, c_dm, n, gnfw..., z=z)
                unitless_test(model, unitless)
            end
        end

        # @test_throws BayesJ.PriorError EinastoModel(mass, fg, c_dm, 2.1, gnfw..., z=z)
        @test_throws BayesJ.PriorError EinastoModel(mass, fg, c_dm, 53 / 3, gnfw..., z=z)
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
            57.0u"kpc",
            3.88,
            1.42u"Mpc",
            0.12,
            5.00,
            10.0,
        ) # Universal values from Arnaud 2010
        unitless_params = Tuple(strip_units(i) for i in params)
        z = 0.5

        model = Vikhlinin2006Model(params..., z=z)
        test_model(model)

        unitless = Vikhlinin2006Model(unitless_params..., z=z)
        unitless_test(model, unitless)
    end

end

function test_piecewise()
    @testset "Piecewise" begin

        params = (
            0.0, 4.0e-23, 5.0,
            10, 3.0e-23, 4.0,
            100, 3e-24, 2.0,
            1000, 0, 0
        )

        z = 0.5

        model = PiecewiseModel(params..., z=z)
        test_model(model)

        @testset "Interpolation" begin
            # TEMPERATURE
            # Values at interpolation points
            @test temperature(model, Quantity(params[1], u"kpc")) == Quantity(params[3], u"keV")
            @test temperature(model, Quantity(params[4], u"kpc")) == Quantity(params[6], u"keV")
            @test temperature(model, Quantity(params[7], u"kpc")) == Quantity(params[9], u"keV")
            @test temperature(model, Quantity(params[10], u"kpc")) == Quantity(params[12], u"keV")
            # Boundary behaviour
            @test temperature(model, Quantity(params[1], u"kpc") - 1u"kpc") == Quantity(params[3], u"keV")
            @test temperature(model, Quantity(params[10], u"kpc") + 1u"kpc") == Quantity(0, u"keV")
            # Interpolation
            @test temperature(model, Quantity(params[1] + params[4], u"kpc") / 2) == Quantity((params[3] + params[6]) / 2, u"keV")

            # DENSITY
            # Values at interpolation points
            @test density(model, Quantity(params[1], u"kpc")) == Quantity(params[2], u"Msun/kpc^3")
            @test density(model, Quantity(params[4], u"kpc")) == Quantity(params[5], u"Msun/kpc^3")
            @test density(model, Quantity(params[7], u"kpc")) == Quantity(params[8], u"Msun/kpc^3")
            @test density(model, Quantity(params[10], u"kpc")) == Quantity(params[11], u"Msun/kpc^3")
            # Boundary behaviour
            @test density(model, Quantity(params[1], u"kpc") - 1u"kpc") == Quantity(params[2], u"Msun/kpc^3")
            @test density(model, Quantity(params[10], u"kpc") + 1u"kpc") == Quantity(0, u"Msun/kpc^3")
            # Interpolation
            @test density(model, Quantity(params[1] + params[4], u"kpc") / 2) == Quantity((params[2] + params[5]) / 2, u"Msun/kpc^3")
        end
    end

end

@testset "Cluster Models" begin
    test_nfw()
    test_einasto()
    test_vikhlinin2006()
    test_piecewise()
end