using BayesJ
using Random, PoissonRandom
using Unitful, UnitfulAstro, DimensionfulAngles
using LinearAlgebra: I

function test_constant_consistency(T::Unitful.Energy=0.5u"keV", ρ::Unitful.Density=1.0e-20u"g/cm^3")
    r::Unitful.Length = 0.5u"kpc"
    z = 0.1
    shape = (3, 3)

    pixel_edge_angle = 0.492u"arcsecondᵃ"
    energy_bins = range(0.7u"keV", 7.0u"keV", length=700)
    n_energy_bins = length(energy_bins) - 1
    exposure_time = 30000.0u"s"
    response_function = 1u"cm^2" * Matrix(I, (n_energy_bins, n_energy_bins))
    centre_radius = 0
    integration_limit = 2 * r

    emission_model = BayesJ.prepare_model_mekal(
        2.0e20u"cm^-2",
        energy_bins,
        z,
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

    function run_sampler(obs, bg, priors, log_dir=nothing)
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
            log_dir=log_dir,
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

    observation, background = make_predicted(r, T, ρ)

    tu = ustrip(u"keV", T)
    ρu = ustrip(u"g/cm^3", ρ)
    label = "$(tu)keV_ρ=$(ρu)gcm3"
    base_log_dir = joinpath("logs", "single_cell3x3_Irmf", label)
    mkpath(base_log_dir)

    # Fit density
    priors = [
        BayesJ.DeltaPrior("x0", 0.0),
        BayesJ.DeltaPrior("y0", 0.0),
        BayesJ.DeltaPrior("r", ustrip(u"Mpc", r)),
        BayesJ.DeltaPrior("T", ustrip(u"keV", T)),
        BayesJ.LogUniformPrior("ρ", 1.0e-29, 1.0e-23),
    ]
    lb, ub, ml, mu = run_sampler(observation, background, priors, joinpath(base_log_dir, "density"))
    if lb[1] < ρu < ub[1]
        BayesJ.@mpiinfo "Density fit passed" lb[1] ρu ub[1]
    else
        BayesJ.@mpiwarn "Density fit failed" lb[1] ρu ub[1]
    end

    # Fit temperature
    priors = [
        BayesJ.DeltaPrior("x0", 0.0),
        BayesJ.DeltaPrior("y0", 0.0),
        BayesJ.DeltaPrior("r", ustrip(u"Mpc", r)),
        BayesJ.UniformPrior("T", 0.0, 8.0),
        BayesJ.DeltaPrior("ρ", ustrip(u"g/cm^3", ρ)),
    ]
    lb, ub, ml, mu = run_sampler(observation, background, priors, joinpath(base_log_dir, "temperature"))
    if lb[1] < tu < ub[1]
        BayesJ.@mpiinfo "Temperature fit passed" lb[1] tu ub[1]
    else
        BayesJ.@mpiwarn "Temperature fit failed" lb[1] tu ub[1]
    end

    # Fit both
    priors = [
        BayesJ.DeltaPrior("x0", 0.0),
        BayesJ.DeltaPrior("y0", 0.0),
        BayesJ.DeltaPrior("r", ustrip(u"Mpc", r)),
        BayesJ.UniformPrior("T", 0.0, 8.0),
        BayesJ.LogUniformPrior("ρ", 1.0e-29, 1.0e-23),
    ]
    lb, ub, ml, mu = run_sampler(observation, background, priors, joinpath(base_log_dir, "temperature_density"))
    if lb[1] < tu < ub[1] && lb[2] < ρu < ub[2]
        BayesJ.@mpiinfo "Temperature and density fit passed" lb[1] tu ub[1] lb[2] ρu ub[2]
    else
        BayesJ.@mpiwarn "Temperature and density fit failed" lb[1] tu ub[1] lb[2] ρu ub[2]
    end
end

for T in [0.5u"keV", 1.0u"keV", 2.0u"keV", 5.0u"keV"]
    for ρ in [1.0e-28u"g/cm^3", 1.0e-27u"g/cm^3", 1.0e-26u"g/cm^3", 1.0e-25u"g/cm^3", 1.0e-24u"g/cm^3"]
        test_constant_consistency(T, ρ)
    end
end

# test_nfw_consistency()