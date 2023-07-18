using Unitful
using Plots

using PyCall
ultranest = pyimport_conda("ultranest", "ultranest", "conda-forge")

export sample

include("gas_models.jl")
include("io.jl")
include("likelihood.jl")

"""
    sample(observed, observed_background, response_function, transform, obs_exposure_time, bg_exposure_time; emission_model, pixel_edge_angle)

Configure some necessary variables and launch ultranest.

* The observed array includes the background.
* The response function includes both the RMF and ARF, as described in `apply_response_function`.
* The emission model should be a function compatible with the requirements of the `surface_brightness` function, which it will be passed to.
* The pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds.
 This area is assumed to be square with the edge angle giving the side length.
* The average effective area is the effective area of the telescope averaged across energies,
 used with the total background rate across all channels (counts per unit telescope area per sky angle per second)
 to calculate the background counts per second per channel per pixel.
"""
function sample(
    observed::T,
    observed_background::T,
    response_function::Matrix,
    transform::Function,
    obs_exposure_time::Unitful.Time,
    bg_exposure_time::Unitful.Time,
    redshift::Real;
    emission_model,
    pixel_edge_angle=0.492u"arcsecond",
    background_rate=8.6e-2u"m^-2/arcminute^2/s",
    average_effective_area=250u"cm^2"
) where {T<:AbstractArray}
    @debug "Preparing for ultranest"

    predicted_bg_rate = background_rate / size(observed)[1] * average_effective_area * pixel_edge_angle^2
    predicted_obs_bg = predicted_bg_rate * obs_exposure_time # Used for adding background to observations
    predicted_bg_bg = predicted_bg_rate * bg_exposure_time # Used for log likelihood

    log_obs_factorial = log_factorial.(observed) + log_factorial.(observed_background)

    @assert all(isfinite, observed)
    @assert all(isfinite, observed_background)
    @assert all(isfinite, log_obs_factorial)

    @assert size(observed) == size(observed_background)

    shape = [i for i in size(observed)][2:3]

    @debug "Observation has shape $(size(observed))"
    @debug "Background has shape $(size(observed_background))"
    @debug "Response matrix has shape $(size(response_function))"

    # a wrapper to handle running the gas model and likelihood calculation
    @debug "Generating likelihood wrapper"
    function likelihood_wrapper(params)
        @debug "Likelihood wrapper called"

        n, _ = size(params)

        params_rows = Vector{Vector{eltype(params)}}(undef, n)
        for i in 1:n
            params_rows[i] = params[i, :]
        end

        predicted = [Model_NFW_GNFW(params_rows[i]...,
            1.062,
            5.4807,
            0.3292,
            1.156,
            redshift,
            shape,
            pixel_edge_angle,
            emission_model,
            obs_exposure_time,
            response_function
        ) .+ predicted_obs_bg for i in 1:n]


        @debug "Predicted results generated"
        # [display(heatmap(dropdims(sum(p, dims=1), dims=1))) for p in predicted]

        return log_likelihood.(
            Ref(observed),
            Ref(observed_background),
            predicted,
            predicted_bg_bg,
            Ref(log_obs_factorial)
        )
    end

    # ultranest setup
    @debug "Creating sampler"
    paramnames = ["MT_200", "fg_200"] # move to pairs with prior objects?
    sampler = ultranest.ReactiveNestedSampler(
        paramnames,
        likelihood_wrapper,
        transform=transform,
        vectorized=true,
        log_dir="logs"
    )

    # run Ultranest
    @info "Launching sampler"
    results = sampler.run()

    # output data
    @debug "Sampler done"
    # print("result has these keys:", keys(results), "\n")
    sampler.print_results()
    sampler.plot_corner()

    return (sampler, results)
end

"""
    sample(data::BayesXDataset, energy_range, priors)

Run Bayesian inference on a given set of `data`, considering only the selected
energy range. An gas emission model `(density, temperature) → emissivity` can be provided.
"""
function sample(
    data::Dataset,
    energy_range,
    priors;
    nHcol=2.2, # units of 10²² atoms per cm⁻²
    redshift=0.1
)
    @info "Loading data"

    observation, observed_background = load_data(data)

    obs = bin_events(observation.first, energy_range, 2000:100:4000, 2000:100:4000)
    bg = bin_events(observed_background.first, energy_range, 2000:100:4000, 2000:100:4000)

    @assert size(obs) == size(bg)

    transform = make_cube_transform(priors...)

    response_function = load_response(data, energy_range)

    @info "Generating emissions model"

    emission_model = prepare_model_mekal(nHcol, 0.1, LinRange(energy_range[1], energy_range[2], size(response_function)[2] + 1)) # we need this +1 but it seems to be one element too short

    sample(obs, bg, response_function, transform, observation.second, observed_background.second, redshift; emission_model=emission_model, pixel_edge_angle=data.pixel_edge_angle)
end