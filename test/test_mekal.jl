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

        T = 1.0u"keV"
        nH = 1.0e-4u"cm^-3"

        direct_call = BayesJ.call_mekal(n_energy_bins, min_energy, max_energy, bin_sizes, ustrip(Cfloat, u"keV", T), ustrip(Cfloat, u"cm^-3", nH))

        # Verify that we return zero emissions in regions of zero gas density
        # or zero temperature.
        @testset "Zero handling" begin
            emission_model = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                0.1,
                use_interpolation=false
            )
            @test all(iszero, emission_model(1.0u"keV", 0.0u"cm^-3"))
            @test all(iszero, emission_model(0.0u"keV", 1.0u"cm^-3"))
            @test all(iszero, emission_model(0.0u"keV", 0.0u"cm^-3"))

            @test all(iszero, BayesJ.call_mekal(energy_bins, 0, 0))
        end

        # Verify that the helper wrapper returns the same results as doing the prep ourselves
        @testset "Wrapper" begin
            @test BayesJ.call_mekal(energy_bins, ustrip(u"keV", T), ustrip(u"cm^-3", nH)) == direct_call
        end

        # Verify that absorption behaves as expected
        @testset "Absorption" begin
            # Verify that zero hydrogen column density means no absorption
            emission_model = BayesJ.prepare_model_mekal(
                0.0u"cm^-2",
                energy_bins,
                0.1,
                use_interpolation=false
            )
            @test_broken emission_model(T, nH) == direct_call

            # Verify that nonzero hydrogen column density means absorption
            emission_model = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                0.1,
                use_interpolation=false
            )
            @test emission_model(T, nH) < direct_call
        end

        # Verify that redshift affects results
        @testset "Redshift" begin
            emission_model_0 = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                0.0,
                use_interpolation=false
            )
            emission_model_01 = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                0.1,
                use_interpolation=false
            )
            emission_model_1 = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                1.0,
                use_interpolation=false
            )
            emission_model_2 = BayesJ.prepare_model_mekal(
                2.2e20u"cm^-2",
                energy_bins,
                2.0,
                use_interpolation=false
            )

            @test emission_model_0(T, nH) != emission_model_01(T, nH)
            @test emission_model_0(T, nH) != emission_model_1(T, nH)
            @test emission_model_0(T, nH) != emission_model_2(T, nH)

            @test emission_model_01(T, nH) != emission_model_1(T, nH)
            @test emission_model_01(T, nH) != emission_model_2(T, nH)

            @test emission_model_1(T, nH) != emission_model_2(T, nH)
        end

        # TODO: Store some results in a data file and check runtime results match them

    end
end

test_mekal()