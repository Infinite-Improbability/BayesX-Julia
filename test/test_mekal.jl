push!(LOAD_PATH, "src/")
using Unitful

# TODO: Test output matches known good output outside of trivial or edge cases

function test_mekal()
    @testset "MEKAL" begin
        # Prepare arguments for call_mekal
        energy_bins = range(0u"keV", 13u"keV", length=10000)
        n_energy_bins = length(energy_bins) - 1
        min_energy = ustrip.(Cfloat, u"keV", energy_bins[1:end-1])
        max_energy = ustrip.(Cfloat, u"keV", energy_bins[2:end])
        abundances = ones(Cfloat, 15)
        nH = BayesJ.HydrogenDensity(abundances)

        T = 1.0u"keV"
        ρ = 1.0e-6u"g/cm^3"

        # Prepare output
        flux = zeros(Cfloat, n_energy_bins)
        direct_call = zeros(Cfloat, n_energy_bins)

        BayesJ.call_mekal(
            direct_call,
            abundances,
            n_energy_bins,
            min_energy,
            max_energy,
            ustrip(Cfloat, u"keV", T),
            ustrip(Cfloat, u"cm^-3", nH(ρ))
        )

        @test all(isfinite, direct_call)
        @test any(!iszero, direct_call)

        # Test that initial flux values don't affect results
        @testset "Nonzero inital flux" begin
            flux .= 3.0
            BayesJ.call_mekal(flux, abundances, energy_bins, ustrip(u"keV", T), ustrip(u"cm^-3", nH(ρ)))
            @test flux == direct_call
        end

        # Verify that the helper wrapper returns the same results as doing the prep ourselves
        @testset "Wrapper" begin
            BayesJ.call_mekal(flux, abundances, energy_bins, ustrip(u"keV", T), ustrip(u"cm^-3", nH(ρ)))
            @test flux == direct_call
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
            fill!(flux, 1.0)
            BayesJ.volume_emissivity!(flux, emission_model, 1.0u"keV", 0.0u"g/cm^3")
            @test all(iszero, flux)
            fill!(flux, 1.0)
            BayesJ.volume_emissivity!(flux, emission_model, 0.0u"keV", 1.0e-23u"g/cm^3")
            @test all(iszero, flux)
            fill!(flux, 1.0)
            BayesJ.volume_emissivity!(flux, emission_model, 0.0u"keV", 0.0u"g/cm^3")
            @test all(iszero, flux)
            fill!(flux, 1.0)
            BayesJ.call_mekal(flux, abundances, energy_bins, 0, 0)
            @test all(iszero, flux)
        end

        # Verify that absorption behaves as expected
        @testset "Absorption" begin
            # Verify that zero hydrogen column density means no absorption
            emission_model = BayesJ.prepare_model_mekal(
                0.0u"cm^-2",
                energy_bins,
                0.0,
            )
            BayesJ.volume_emissivity!(flux, emission_model, T, ρ)
            @test_broken flux == direct_call # this no longer succeeds because the application of mekal coefficents was moved out of call_mekal and into volume_emissivity!

            # Verify that nonzero hydrogen column density means absorption
            emission_model = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                0.0,
            )
            BayesJ.volume_emissivity!(flux, emission_model, T, ρ)
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

            BayesJ.volume_emissivity!(flux0, emission_model_0, T, ρ)
            BayesJ.volume_emissivity!(flux01, emission_model_01, T, ρ)
            BayesJ.volume_emissivity!(flux1, emission_model_1, T, ρ)
            BayesJ.volume_emissivity!(flux2, emission_model_2, T, ρ)

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


BayesJ.download_spectral_fitting_data(; progress=true, verbose=true)
test_mekal()