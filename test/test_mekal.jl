push!(LOAD_PATH, "src/")
using Unitful

function test_mekal()
    @testset "MEKAL" begin
        # Prepare arguments for call_mekal
        energy_bins = range(0u"keV", 13u"keV", length=10000)
        n_energy_bins = length(energy_bins) - 1
        min_energy = ustrip.(Cfloat, u"keV", energy_bins[1:end-1])
        max_energy = ustrip.(Cfloat, u"keV", energy_bins[2:end])
        bin_sizes = max_energy - min_energy
        abundances = ones(Cfloat, 15)

        T = 1.0u"keV"
        nH = 1.0e-4u"cm^-3"

        # Prepare output
        flux = zeros(Cfloat, n_energy_bins)
        direct_call = zeros(Cfloat, n_energy_bins)

        BayesJ.call_mekal(
            direct_call,
            abundances,
            n_energy_bins,
            min_energy,
            max_energy,
            bin_sizes,
            ustrip(Cfloat, u"keV", T),
            ustrip(Cfloat, u"cm^-3", nH)
        )

        @test all(isfinite, direct_call)

        # Test that initial flux values don't affect results
        @testset "Nonzero inital flux" begin
            flux .= 3.0
            BayesJ.call_mekal(flux, abundances, energy_bins, ustrip(u"keV", T), ustrip(u"cm^-3", nH))
            @test flux == direct_call
        end

        # Verify that the helper wrapper returns the same results as doing the prep ourselves
        @testset "Wrapper" begin
            BayesJ.call_mekal(flux, abundances, energy_bins, ustrip(u"keV", T), ustrip(u"cm^-3", nH))
            @test flux == direct_call
            display(flux ≈ direct_call)
        end

        # Verify that we return zero emissions in regions of zero gas density
        # or zero temperature.
        @testset "Zero handling" begin
            emission_model = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                0.1,
                abundances
            )
            emission_model(flux, 1.0u"keV", 0.0u"cm^-3")
            @test all(iszero, flux)
            emission_model(flux, 0.0u"keV", 1.0u"cm^-3")
            @test all(iszero, flux)
            emission_model(flux, 0.0u"keV", 0.0u"cm^-3")
            @test all(iszero, flux)

            BayesJ.call_mekal(flux, abundances, energy_bins, 0, 0)
            @test all(iszero, flux)
        end

        # Verify that absorption behaves as expected
        @testset "Absorption" begin
            # Verify that zero hydrogen column density means no absorption
            emission_model = BayesJ.prepare_model_mekal(
                0.0u"cm^-2",
                energy_bins,
                0.1,
            )
            emission_model(flux, T, nH)
            @test_broken flux == direct_call

            # Verify that nonzero hydrogen column density means absorption
            emission_model = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                0.1,
            )
            emission_model(flux, T, nH)
            @test flux < direct_call
        end

        # Verify that redshift affects results
        @testset "Redshift" begin
            emission_model_0 = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                0.0,
                abundances
            )
            emission_model_01 = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                0.1,
                abundances
            )
            emission_model_1 = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                1.0,
                abundances
            )
            emission_model_2 = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                2.0,
                abundances
            )

            flux0 = copy(flux)
            flux01 = copy(flux)
            flux1 = copy(flux)
            flux2 = copy(flux)

            emission_model_0(flux0, T, nH)
            emission_model_01(flux01, T, nH)
            emission_model_1(flux1, T, nH)
            emission_model_2(flux2, T, nH)

            @test flux0 != flux01
            @test flux0 != flux1
            @test flux0 != flux2

            @test flux01 != flux1
            @test flux01 != flux2

            @test flux1 != flux2
        end

        # TODO: Store some results in a data file and check runtime results match them

    end
end

test_mekal()