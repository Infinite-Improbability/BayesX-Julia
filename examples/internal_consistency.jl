using BayesJ
using Random, PoissonRandom
using Unitful, UnitfulAstro, DimensionfulAngles
using LinearAlgebra: I
using CairoMakie

function plot(obs::Array, pred::Array, path::AbstractString)
    if BayesJ.isroot()
        if pred[:, 1, 1] != pred[:, 4, 4]
            @warn "Something funky with the best fit array"
        end
        slice_obs = Vector{Float64}(obs[:, 4, 4])
        slice_pred = Vector{Float64}(pred[:, 4, 4])

        f = Figure()
        ax = Axis(f[1, 1], xlabel="Channel", ylabel="Counts", yscale=Makie.pseudolog10)
        lines!(slice_obs ./ maximum(slice_obs), label="Observation")
        lines!(slice_pred ./ maximum(slice_pred), label="Best Fit")
        axislegend()

        save(joinpath(path, "spectra.svg"), f)
    end
end

function test_constant_consistency(T::Unitful.Energy, ρ::Unitful.Density)
    tu = ustrip(u"keV", T)
    ρu = ustrip(u"g/cm^3", ρ)
    label = "$(tu)keV_ρ=$(ρu)gcm3"
    base_log_dir = joinpath("..", "logs", "constant", "acisi-cy0", label)

    if isfile(joinpath(base_log_dir, "temperature_density", "run1", "info", "results.json"))
        BayesJ.@mpiwarn "Skipping run as folder exists for this combination" T ρ base_log_dir
        return
    end

    mkpath(base_log_dir)

    r::Unitful.Length = 5.0u"Mpc"
    z = 0.1
    shape = (9, 9)

    pixel_edge_angle = 20.0u"arcsecondᵃ"
    energy_bins = range(0.7u"keV", 7.0u"keV", step=0.05u"keV")
    exposure_time = 3.0e6u"s"

    data = FITSData(
        "",
        "",
        "../data/tng/acisi_aimpt_cy0.arf",
        "../data/tng/acisi_aimpt_cy0.rmf",
        0.492u"arcsecondᵃ"
    )

    response_function, energy_bins, _ = BayesJ.load_response(data, 0.7u"keV", 7.0u"keV")
    centre_radius = 0
    integration_limit = 2 * r

    emission_model = BayesJ.prepare_model_mekal(
        0.0e0u"cm^-2",
        energy_bins,
        z,
    )

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

        try
            output_dir = sampler.logs["run_dir"]
            if output_dir isa AbstractString
                plot(obs, best_fit_observation, joinpath(output_dir, "plots"))
            end
        catch e
            if !(e isa KeyError)
                rethrow(e)
            end
            BayesJ.@mpidebug "Skipping spectra plot because unable to find log information"
        end

        lower_bound = result["posterior"]["errlo"]
        upper_bound = result["posterior"]["errup"]
        mean = result["posterior"]["mean"]
        std = result["posterior"]["stdev"]

        mean_lower = mean - 2 * std
        mean_upper = mean + 2 * std

        return lower_bound, upper_bound, mean_lower, mean_upper
    end

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
    bg_rate = rand() * min(1, maximum(predicted_count_rate))
    background = Array{Int64}(undef, size(observation)...)
    for i in eachindex(background)
        background[i] = pois_rand(bg_rate) + 1
    end

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
    for ρ in [3.0e-28u"g/cm^3", 3.0e-27u"g/cm^3", 3.0e-26u"g/cm^3"]
        test_constant_consistency(T, ρ)
        # GC.gc(true)
        # ccall(:malloc_trim, Int32, (Int32,), 0)
    end
end