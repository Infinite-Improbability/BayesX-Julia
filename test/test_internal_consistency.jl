using PoissonRandom
using PyCall
using Random
using Test
using Unitful, UnitfulAstro, DimensionfulAngles
using BayesJ

numpy = pyimport("numpy")

numpy.random.seed(42)
Random.seed!(4242)

function test_single_cell_consistency()
    numpy.random.seed(42)
    Random.seed!(4242)

    r = 0.5u"kpc"
    T = 0.5u"keV"
    ρ = 1.0e-20u"g/cm^3"

    z = 0.1
    shape = (1, 1)

    pixel_edge_angle = 0.0492u"arcsecondᵃ"
    energy_bins = range(0.7u"keV", 7.0u"keV", length=200)
    exposure_time = 3.0u"s"
    response_function = rand(Float64, (100, length(energy_bins) - 1)) * 1u"cm^2"
    centre_radius = 0
    integration_limit = 1.0u"kpc"

    emission_model = BayesJ.prepare_model_mekal(
        2.0e20u"cm^-2",
        energy_bins,
        z,
        use_interpolation=false
    )

    function make_predicted(r, T, ρ)
        temperature, density = Model_NFW(r, T, ρ)

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

        mean_lower = mean - 2 * std
        mean_upper = mean + 2 * std

        return lower_bound, upper_bound, mean_lower, mean_upper
    end

    @testset "Single Cell IC" begin

        @testset "Single Cell IC (fit ρ)" begin
            priors = [
                BayesJ.DeltaPrior("x0", 0.0),
                BayesJ.DeltaPrior("y0", 0.0),
                BayesJ.DeltaPrior("r", ustrip(u"Mpc", r)),
                BayesJ.DeltaPrior("T", ustrip(u"keV", T)),
                BayesJ.LogUniformPrior("ρ", 1.0e-24, 1.0e-15),
            ]
            observation, background = make_predicted(r, T, ρ)
            lower_bound, upper_bound, mean_lower, mean_upper = run_sampler(observation, background, priors)
            @test lower_bound[1] < ustrip(u"g/cm^3", ρ) < upper_bound[1] || mean_lower[1] < ustrip(u"g/cm^3", ρ) < mean_upper[1]
        end

        @testset "Single Cell IC (fit T)" begin
            priors = [
                BayesJ.DeltaPrior("x0", 0.0),
                BayesJ.DeltaPrior("y0", 0.0),
                BayesJ.DeltaPrior("r", ustrip(u"Mpc", r)),
                BayesJ.UniformPrior("T", 0.0, 10.0),
                BayesJ.DeltaPrior("ρ", ustrip(u"g/cm^3", ρ)),
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

function test_nfw_consistency()

    mass = 3e14
    fg = 0.13
    c_dm = 3.0
    gnfw = [1.0510, 5.4905, 0.3081, 1.177] # Universal values from Arnaud 2010

    z = 0.1
    shape = (50, 50)

    pixel_edge_angle = 0.0492u"arcsecondᵃ"
    energy_bins = range(0.7u"keV", 7.0u"keV", length=200)
    exposure_time = 3.0e5u"s"
    response_function = rand(Float64, (100, length(energy_bins) - 1)) * 1u"cm^2"
    centre_radius = 0
    integration_limit = 10u"Mpc"

    emission_model = BayesJ.prepare_model_mekal(
        2.0e20u"cm^-2",
        energy_bins,
        z,
        use_interpolation=false
    )

    temperature, density = Model_NFW(mass, fg, c_dm, gnfw..., z=z)

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
            cluster_model=Model_NFW,
            emission_model=emission_model,
            param_wrapper=param_wrapper,
            pixel_edge_angle=pixel_edge_angle,
            centre_radius=centre_radius,
            log_dir=nothing,
            integration_limit=integration_limit,
            ultranest_run_args=(
                max_num_improvement_loops=3,
                min_num_live_points=100,
                # show_status=false,
                # viz_callback=false
            )
        )

        lower_bound = result["posterior"]["errlo"]
        upper_bound = result["posterior"]["errup"]
        mean = result["posterior"]["mean"]
        std = result["posterior"]["stdev"]

        mean_lower = mean - 2 * std
        mean_upper = mean + 2 * std

        return lower_bound, upper_bound, mean_lower, mean_upper
    end

    @testset "Single Cell IC" begin

        @testset "NFW IC (fit Mt, fg)" begin
            priors = [
                DeltaPrior("x0", 0.0),
                DeltaPrior("y0", 0.0),
                UniformPrior("M_T", 1e14, 1e15),
                UniformPrior("f_g", 0.0, 0.5),
                DeltaPrior("c_500", c_dm),
                DeltaPrior("a", gnfw[1]),
                DeltaPrior("b", gnfw[2]),
                DeltaPrior("c", gnfw[3]),
                DeltaPrior("c_500_GNFW", gnfw[4]
                )
            ]
            lower_bound, upper_bound, mean_lower, mean_upper = run_sampler(observation, background, priors)
            @test lower_bound[1] < mass < upper_bound[1] || mean_lower[1] < mass < mean_upper[1]
            @test lower_bound[2] < fg < upper_bound[2] || mean_lower[2] < fg < mean_upper[2]
        end
    end

end

test_single_cell_consistency()
# test_nfw_consistency()