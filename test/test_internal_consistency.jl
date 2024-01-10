using PoissonRandom

# include("../src/BayesJ.jl") # for completion during dev


function test_single_cell_consistency()
    r = 0.5u"kpc"
    T = 2.0u"keV"
    ρ = 1.0e-20u"g/cm^3"

    z = 0.1
    shape = (1, 1)

    pixel_edge_angle = 0.0492u"arcsecondᵃ"
    energy_bins = range(0.7u"keV", 3.0u"keV", length=100)
    exposure_time = 3.0u"s"
    response_function = rand(Float64, (50, length(energy_bins) - 1)) * 1u"cm^2"
    centre_radius = 0
    integration_limit = 1.0u"kpc"

    emission_model = BayesJ.prepare_model_mekal(
        2.0e20u"cm^-2",
        energy_bins,
        z,
        use_interpolation=false
    )

    function make_predicted(r, T, ρ)
        temperature, density = Model_Constant(r, T, ρ)

        predicted_count_rate = BayesJ.make_observation(
            temperature,
            density,
            z,
            shape,
            pixel_edge_angle,
            emission_model,
            exposure_time,
            response_function,
            (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
            centre_radius,
            limit=integration_limit,
        )

        @assert all(isfinite, predicted_count_rate)

        observation = pois_rand.(predicted_count_rate)
        background_rate = rand()
        background = Array{Int64}(undef, size(observation)...)
        for i in eachindex(background)
            background[i] = pois_rand(background_rate) + 1
        end

        return observation, background
    end

    function run_sampler(obs, bg, priors)
        prior_transform, param_wrapper = BayesJ.make_cube_transform(priors...)
        prior_names = [p.name for p in priors if !isa(p, DeltaPrior)]

        sampler, result, best_fit_observation = BayesJ.sample(
            obs + bg,
            bg,
            response_function,
            prior_transform,
            exposure_time,
            exposure_time, # TODO: Make different
            z;
            prior_names=prior_names,
            cluster_model=Model_Constant,
            emission_model=emission_model,
            param_wrapper=param_wrapper,
            pixel_edge_angle=pixel_edge_angle,
            centre_radius=centre_radius,
            log_dir=nothing,
            integration_limit=integration_limit,
            ultranest_run_args=(
                max_num_improvement_loops=3,
                min_num_live_points=100,
                show_status=false,
                viz_callback=false
            )
        )

        lower_bound = result["posterior"]["errlo"]
        upper_bound = result["posterior"]["errup"]
        mean = result["posterior"]["mean"]
        std = result["posterior"]["stdev"]

        mean_lower = mean - std
        mean_upper = mean + std

        return lower_bound, upper_bound, mean_lower, mean_upper
    end

    @testset "Single Cell IC (fit ρ)" begin
        priors = [
            BayesJ.DeltaPrior("x0", 0.0),
            BayesJ.DeltaPrior("y0", 0.0),
            BayesJ.DeltaPrior("r", ustrip(u"Mpc", r)),
            # BayesJ.UniformPrior("T", 0.0, 10.0),
            BayesJ.DeltaPrior("T", 2.0),
            BayesJ.LogUniformPrior("ρ", 1.0e-24, 1.0e-15),
        ]
        observation, background = make_predicted(r, T, ρ)
        lower_bound, upper_bound, mean_lower, mean_upper = run_sampler(observation, background, priors)
        @test lower_bound[1] < ustrip(u"g/cm^3", ρ) < upper_bound[1] || mean_lower[1] < ustrip(u"g/cm^3", ρ) < mean_upper[1]
    end

    @testset "Single Cell IC" begin

        @testset "Single Cell IC (fit T)" begin
            priors = [
                BayesJ.DeltaPrior("x0", 0.0),
                BayesJ.DeltaPrior("y0", 0.0),
                BayesJ.DeltaPrior("r", ustrip(u"Mpc", r)),
                BayesJ.UniformPrior("T", 0.0, 10.0),
                BayesJ.DeltaPrior("ρ", 1.0e-20),
            ]
            observation, background = make_predicted(r, T, ρ)
            lower_bound, upper_bound, mean_lower, mean_upper = run_sampler(observation, background, priors)
            @test lower_bound[1] < ustrip(u"keV", T) < upper_bound[1] || mean_lower[1] < ustrip(u"keV", T) < mean_upper[1]
        end

        @testset "Single Cell IC (fit T and ρ)" begin
            priors = [
                BayesJ.DeltaPrior("x0", 0.0),
                BayesJ.DeltaPrior("y0", 0.0),
                BayesJ.DeltaPrior("r", ustrip(u"Mpc", r)),
                BayesJ.UniformPrior("T", 0.0, 10.0),
                BayesJ.LogUniformPrior("ρ", 1.0e-24, 1.0e-15),
            ]
            observation, background = make_predicted(r, T, ρ)
            lower_bound, upper_bound, mean_lower, mean_upper = run_sampler(observation, background, priors)
            @test lower_bound[1] < ustrip(u"keV", T) < upper_bound[1] || mean_lower[1] < ustrip(u"keV", T) < mean_upper[1]
            @test lower_bound[2] < ustrip(u"g/cm^3", ρ) < upper_bound[2] || mean_lower[2] < ustrip(u"g/cm^3", ρ) < mean_upper[2]
        end
    end

end

test_single_cell_consistency()