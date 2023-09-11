using Unitful, DimensionfulAngles

using PyCall
ultranest = pyimport_conda("ultranest", "ultranest", "conda-forge")
stepsampler = pyimport_conda("ultranest.stepsampler", "ultranest")

export sample

include("mpi.jl")
include("cluster_models.jl")
include("io.jl")
include("likelihood.jl")

"""
    sample(observed, observed_background, response_function, transform, obs_exposure_time, bg_exposure_time, redshift; prior_names, cluster_model, emission_model, param_wrapper, pixel_edge_angle, background_rate, average_effective_area)

Configure some necessary variables and launch ultranest.

* The observed array includes the background.
* The response function includes both the RMF and ARF, as described in `apply_response_function`.
* The emission model should be a function compatible with the requirements of the `surface_brightness` function, which it will be passed to.
* The pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds.
 This area is assumed to be square with the edge angle giving the side length.
* The average effective area is the effective area of the telescope averaged across energies,
 used with the total background rate across all channels (counts per unit telescope area per sky angle per second)
 to calculate the background counts per second per channel per pixel.
* The first two priors should always be "x0" and "y0", giving cluster centre position and the order of prior_names must match the transform.
* The cluster model should be a function that takes the parameters (and redshift as a kwarg) and returns `(gas_temperature,gas_density)` as functions of
radius which return their respective quantities with units.
* `param_wrapper` takes the output of the transform function and adds any additional arguements necessary for the model.
"""
function sample(
    observed::T,
    observed_background::T,
    response_function::Matrix,
    transform::Function,
    obs_exposure_time::Unitful.Time,
    bg_exposure_time::Unitful.Time,
    redshift::Real;
    prior_names::Vector{<:AbstractString},
    cluster_model::Function,
    emission_model,
    param_wrapper::Function,
    pixel_edge_angle=0.492u"arcsecondᵃ",
    background_rate=8.4e-6u"cm^-2/arcminuteᵃ^2/s",
    average_effective_area=250u"cm^2",
    centre_radius=0,
    mask=nothing,
    use_stepsampler=false,
    log_dir="logs",
    resume="subfolder"
) where {T<:AbstractArray}
    @mpidebug "Preparing for ultranest"

    predicted_bg_rate = background_rate / size(observed)[1] * average_effective_area * pixel_edge_angle^2
    predicted_obs_bg = predicted_bg_rate * obs_exposure_time # Used for adding background to observations
    predicted_bg_bg = predicted_bg_rate * bg_exposure_time # Used for log likelihood

    log_obs_factorial = log_factorial.(observed) + log_factorial.(observed_background)

    @assert predicted_obs_bg > 0
    @assert all(isfinite, observed)
    @assert all(isfinite, observed_background)
    @assert all(isfinite, log_obs_factorial)

    @assert size(observed) == size(observed_background)

    shape = [i for i in size(observed)][2:3]

    @mpiinfo "Observation has shape $(size(observed))"
    @mpiinfo "Background has shape $(size(observed_background))"
    @mpiinfo "Response matrix has shape $(size(response_function))"

    # a wrapper to handle running the gas model and likelihood calculation
    @mpidebug "Generating likelihood wrapper"
    function likelihood_wrapper(params)
        @mpirankeddebug "Likelihood wrapper called" params

        full_params = param_wrapper(params)
        @mpirankeddebug "Full parameters are" full_params

        try
            gas_temperature, gas_density = cluster_model(
                full_params[3:end]...;
                z=redshift
            )

            predicted = make_observation(
                gas_temperature,
                gas_density,
                redshift,
                shape,
                pixel_edge_angle,
                emission_model,
                obs_exposure_time,
                response_function,
                (full_params[1], full_params[2]),
                centre_radius,
                mask=mask
            )

            predicted .= predicted .+ predicted_obs_bg

            @mpirankeddebug "Predicted results generated"

            return log_likelihood(
                observed,
                observed_background,
                predicted,
                predicted_bg_bg,
                log_obs_factorial
            )
        catch e
            if e isa PriorError
                return e.likelihood
            end
            rethrow(e)
        end
    end

    # ultranest setup
    @mpidebug "Creating sampler"
    sampler = ultranest.ReactiveNestedSampler(
        prior_names,
        likelihood_wrapper,
        transform=transform,
        vectorized=false,
        log_dir=log_dir,
        resume=resume
    )

    if use_stepsampler
        @mpidebug "Creating stepsampler"
        sampler.stepsampler = stepsampler.SliceSampler(
            nsteps=2 * length(prior_names),
            generate_direction=stepsampler.generate_mixture_random_direction,
            adaptive_nsteps="move-distance",
            max_nsteps=24,
            region_filter=true
        )
    end

    # run Ultranest
    @mpiinfo "Launching sampler"
    results = sampler.run(
    # region_class=ultranest.mlfriends.RobustEllipsoidRegion
    )

    # output data
    @mpidebug "Sampler done"
    # print("result has these keys:", keys(results), "\n")
    sampler.print_results()
    sampler.plot_corner()
    sampler.plot_trace()
    sampler.plot_run()

    return (sampler, results)
end

"""
    sample(data::Dataset, energy_range::AbstractRange{Unitful.Energy}, cluster_model::Function, priors::AbstractVector{Prior}, nhCol::SurfaceDensity, redshift, x, y)

Run Bayesian inference on a given set of `data` considering only the selected
energy range.

* An gas emission model `(density, temperature) → emissivity` can be provided.
* The first two priors should always be `x0` and `y0`, giving cluster centre position.
* `x` and `y` are tuples of `(min, max)`.
"""
function sample(
    data::Dataset,
    energy_range::AbstractRange{<:Unitful.Energy},
    cluster_model::Function,
    priors::AbstractVector{<:Prior},
    nHcol::SurfaceDensity,
    redshift::Real,
    x::NTuple{2,<:Real},
    y::NTuple{2,<:Real};
    bin_size::Real=10,
    use_interpolation::Bool=true,
    centre_radius=0,
    mask=nothing,
    kwargs...
)
    @argcheck [p.name for p in priors[1:2]] == ["x0", "y0"]

    @mpiinfo "Loading data"
    observation, observed_background = load_data(data)

    x_edges = x[1]:bin_size:x[2]
    y_edges = y[1]:bin_size:y[2]

    @mpidebug "Spatial bounds" x_edges y_edges

    obs = bin_events(data, observation.first, energy_range, x_edges, y_edges)
    bg = bin_events(data, observed_background.first, energy_range, x_edges, y_edges)
    pixel_edge_angle = bin_size * data.pixel_edge_angle
    @assert size(obs) == size(bg)

    @mpidebug "Making transform"
    prior_names = [p.name for p in priors if !isa(p, DeltaPrior)]
    transform, param_wrapper = make_cube_transform(priors...)

    @mpidebug "Calling load_response"
    response_function = load_response(data, energy_range)
    @mpidebug "Response function loaded. Size[2] should be one less than energy range" size(response_function) length(energy_range)
    @assert size(response_function, 2) == (length(energy_range) - 1)

    @mpiinfo "Generating emissions model"
    emission_model = prepare_model_mekal(nHcol, energy_range, redshift, use_interpolation=use_interpolation)

    if mask isa AbstractString
        @mpidebug "Loading mask"
        mask = load_mask(mask, x_edges, y_edges)
        @mpiinfo "Mask size is $(size(mask))"
    end

    sample(
        obs,
        bg,
        response_function,
        transform,
        observation.second,
        observed_background.second,
        redshift;
        prior_names=prior_names,
        cluster_model=cluster_model,
        emission_model=emission_model,
        param_wrapper=param_wrapper,
        pixel_edge_angle=pixel_edge_angle,
        centre_radius=centre_radius,
        mask=mask,
        kwargs...
    )
end