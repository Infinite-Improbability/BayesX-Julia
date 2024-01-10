using PoissonRandom

# include("../src/BayesJ.jl") # for completion during dev

function test_single_cell_consistency()
    T = 2.0u"keV"
    ρ = 1.0e-4u"g/cm^3"
    z = 0.1
    shape = (1, 1)
    pixel_edge_angle = 0.492u"arcsecondᵃ"
    energy_bins = range(0.7u"keV", 3.0u"keV", length=100)
    exposure_time = 3.0e5u"s"
    response_function = rand(Float64, (10, length(energy_bins))) * 1u"cm^2"
    centre_radius = 0

    temperature, density = Model_Constant(T, ρ)
    emission_model = BayesJ.prepare_model_mekal(
        2.0e20u"cm^-2",
        energy_bins,
        z,
        use_interpolation=false
    )

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
        centre_radius
    )

    observation = pois_rand.(predicted_count_rate)
    background_rate = rand()
    background = Matrix{Int64}(undef, size(observation)...)
    for i in eachindex(background)
        background[i] = pois_rand(background_rate)
    end

    priors = [
        BayesJ.DeltaPrior("x0", 0.0),
        BayesJ.DeltaPrior("x0", 0.0),
        BayesJ.UniformPrior(0.0, 10.0),
        BayesJ.UniformPrior(0.0, 10.0)
    ]
    prior_transform, param_wrapper = BayesJ.make_cube_transform(priors...)

    sampler, result, best_fit_observation = BayesJ.sample(
        observation + background,
        background,
        response_function,
        prior_transform,
        exposure_time,
        exposure_time, # TODO: Make different
        z,
        ["T", "ρ"],
        Model_Constant,
        emission_model,
        param_wrapper,
        pixel_edge_angle,
        centre_radius,
        log_dir=nothing,
        ultranest_run_args=(max_num_improvement_loops=3, min_num_live_points=100)
    )

    lower_bound = result["posterior"]["errlo"]
    upper_bound = result["posterior"]["errhi"]

    @testset "Single Cell IC" begin
        @test lower_bound[0] < T < upper_bound[0]
        @test lower_bound[1] < ρ < upper_bound[1]
    end

end

test_single_cell_consistency()