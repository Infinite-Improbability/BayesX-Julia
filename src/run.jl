using Unitful, DimensionfulAngles
using ArgCheck

using PyCall

const ultranest = PyNULL()
const stepsampler = PyNULL()
const mlfriends = PyNULL()

function __init__()
    copy!(ultranest, pyimport_conda("ultranest", "ultranest", "conda-forge"))
    copy!(stepsampler, pyimport_conda("ultranest.stepsampler", "ultranest"))
    copy!(mlfriends, pyimport_conda("ultranest.mlfriends", "ultranest"))
end

export sample, mlfriends

include("mpi.jl")
include("cluster_models.jl")
include("io.jl")
include("likelihood.jl")
include("blobs.jl")


function predict_counts_with_params(
    full_params::P;
    cluster_model::Function,
    emission_model::Function,
    redshift::Real,
    predicted_bg_over_obs_time::B,
    shape::NTuple{2,<:Integer},
    pixel_edge_angle::DimensionfulAngles.Angle,
    observation_exposure_time::Unitful.Time,
    response_function::AbstractMatrix,
    centre_radius::Integer,
    mask::M,
    integration_limit::Unitful.Length
)::Array{Union{Float64,Missing},3} where {P<:Union{AbstractVector,Tuple},B<:Union{AbstractArray,Real},M<:Union{Nothing,AbstractMatrix{Bool}}}
    centre = full_params[1:2]
    model_parameters = full_params[3:end]

    gas_temperature, gas_density = cluster_model(
        model_parameters...;
        z=redshift
    )

    predicted = make_observation(
        gas_temperature,
        gas_density,
        redshift,
        shape,
        pixel_edge_angle,
        emission_model,
        observation_exposure_time,
        response_function,
        centre,
        centre_radius,
        mask=mask,
        limit=integration_limit
    )

    # this intrinsically broadcasts along the energy axis
    predicted .+ predicted_bg_over_obs_time
end

function prepare_background(
    observed::AbstractArray{<:Integer,3},
    observed_background::AbstractArray{<:Integer,3},
    obs_exposure_time::Unitful.Time,
    bg_exposure_time::Unitful.Time,
)::NTuple{2,Vector{Float64}}
    # implicitly includes average effective area and pixel edge angle
    bg_count_rate = [mean(@view observed_background[i, :, :]) for i in axes(observed_background, 1)] ./ bg_exposure_time
    number_of_zeros_in_bg = count(i -> i == 0u"s^-1", bg_count_rate)
    if number_of_zeros_in_bg > 0
        @mpiwarn "Background count has areas of zero counts, replacing with fallback value" number_of_zeros_in_bg
        replace!(bg_count_rate, 0u"s^-1" => 0.1e-6 / bg_exposure_time)
    end

    @assert all(i -> i > 0u"s^-1", bg_count_rate)
    @mpidebug "Background rate estimated" bg_count_rate

    # vector of background as a function of energy
    predicted_obs_bg = bg_count_rate * obs_exposure_time # Used for adding background to observations
    predicted_bg_bg = bg_count_rate * bg_exposure_time # Used for log likelihood

    @assert all(isfinite, predicted_obs_bg)
    @assert all(isfinite, predicted_bg_bg)
    @assert all(i -> i > 0, predicted_obs_bg)
    @assert all(i -> i > 0, predicted_bg_bg)
    @assert length(predicted_bg_bg) == size(observed, 1)

    return predicted_obs_bg, predicted_bg_bg
end

function plateau_test(
    cube_transform::Function,
    prior_names::Vector{<:AbstractString},
    likelihood_wrapper::Function;
    n_tests::Integer=500
)
    @mpiinfo "Testing likelihood calculation"
    test_runs = Matrix{Float64}(undef, length(prior_names) + 1, n_tests)

    @showprogress for i in axes(test_runs, 2)
        params = cube_transform(rand(length(prior_names)))
        test_runs[:, i] .= [params; likelihood_wrapper(params)]
    end
    likelihoods = test_runs[end, :]

    @mpiinfo "Test runs results" n_tests minimum(likelihoods) maximum(likelihoods) mean(likelihoods) median(likelihoods) mode(likelihoods)
    value = mode(likelihoods)
    occurrences = count(i -> i == value, likelihoods)
    if occurrences / n_tests > 0.1
        @mpiwarn "Likelihood plateau detected" value occurrences
    else
        @mpiinfo "Most common likelihood" value occurrences
    end
    plateau_points = test_runs[:, likelihoods.==value]
    other_points = test_runs[:, likelihoods.!==value]
    @assert size(plateau_points, 2) == occurrences
    @assert size(other_points, 2) == n_tests - occurrences
    for i in axes(plateau_points, 1)
        if i == axes(plateau_points, 1)[end]
            continue
        end
        @mpiinfo "Plateau parameter" i minimum(plateau_points[i, :]) maximum(plateau_points[i, :]) mean(plateau_points[i, :]) median(plateau_points[i, :]) mode(plateau_points[i, :])
        @mpiinfo "Other parameter" i minimum(other_points[i, :]) maximum(other_points[i, :]) mean(other_points[i, :]) median(other_points[i, :]) mode(other_points[i, :])
        if minimum(plateau_points[i, :]) > maximum(other_points[i, :]) || maximum(plateau_points[i, :]) < minimum(other_points[i, :])
            @mpiwarn "Plateau and other parameter ranges do not overlap" i
        elseif !any(minimum(plateau_points[i, :]) .< other_points[i, :] .< maximum(plateau_points[i, :]))
            @mpiwarn "Plateau and other parameter ranges do not overlap but plateau is embedded in other" i
        elseif maximum(plateau_points[i, :]) > maximum(other_points[i, :])
            pct = (maximum(plateau_points[i, :]) - maximum(other_points[i, :])) / maximum(other_points[i, :]) * 100
            if pct > 5
                @mpiwarn "Plateau range extends above other range by $pct%."
            end
        elseif minimum(plateau_points[i, :]) < minimum(other_points[i, :])
            pct = (minimum(other_points[i, :]) - minimum(plateau_points[i, :])) / minimum(other_points[i, :]) * 100
            if pct > 5
                @mpiwarn "Plateau range extends below other range by $pct%."
            end
        end
    end
end

"""
    sample(observed,observed_background, response_function, transform, obs_exposure_time, bg_exposure_time, redshift; prior_names, cluster_model, emission_model, param_wrapper, pixel_edge_angle)

Configure some necessary variables and launch ultranest.

* The observed array includes the background. The first dimension is energy, the other two are spatial.
* The response function includes both the RMF and ARF, as described in `apply_response_function`.
* The emission model should be a function compatible with the requirements of the `surface_brightness` function, which it will be passed to.
* The pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds.
 This area is assumed to be square with the edge angle giving the side length.
* The first two priors should always be "x0" and "y0", giving cluster centre position and the order of prior_names must match the transform.
* The cluster model should be a function that takes the parameters (and redshift as a kwarg) and returns `(gas_temperature,gas_density)` as functions of
radius which return their respective quantities with units.
* `param_wrapper` takes the output of the transform function and adds any additional arguements necessary for the model.
* `centre_radius`, `mask` and `integration_limit` are passed through to [`make_observation`](@ref)
 
"""
function sample(
    observed::T,
    observed_background::T,
    response_function::AbstractMatrix,
    transform::Function,
    obs_exposure_time::Unitful.Time,
    bg_exposure_time::Unitful.Time,
    redshift::Real;
    prior_names::AbstractVector{<:AbstractString},
    cluster_model::Function,
    emission_model::Function,
    param_wrapper::Function,
    pixel_edge_angle::DimensionfulAngles.Angle=0.492u"arcsecondᵃ",
    centre_radius::Integer=0,
    mask::Union{Matrix{Bool},Nothing}=nothing,
    integration_limit::Unitful.Length=Quantity(10, u"Mpc"),
    use_stepsampler=false,
    log_dir="logs",
    resume="subfolder",
    plateau_tests=0,
    ultranest_run_args=NamedTuple()
) where {T<:AbstractArray}
    @mpidebug "Preparing for ultranest"

    @argcheck all(isfinite, observed)
    @argcheck all(isfinite, observed_background)
    @argcheck all(i -> i >= 0, observed)
    @argcheck all(i -> i >= 0, observed_background)
    @argcheck size(observed) == size(observed_background)
    @argcheck size(observed, 1) == size(response_function, 1)

    @mpiinfo "Preparing background"
    predicted_obs_bg, predicted_bg_bg = prepare_background(
        observed,
        observed_background,
        obs_exposure_time,
        bg_exposure_time
    )

    log_obs_factorial = log_factorial.(observed) + log_factorial.(observed_background)
    @assert all(isfinite, log_obs_factorial)

    shape = (size(observed, 2), size(observed, 3))

    @mpiinfo "Observation has shape $(size(observed))"
    @mpiinfo "Background has shape $(size(observed_background))"
    @mpiinfo "Response matrix has shape $(size(response_function))"

    if isnothing(mask)
        mask = zeros(Bool, shape)
    end

    # generate count rates matrix for given parameters
    function predict_counts(params::AbstractVector{Float64})
        full_params = param_wrapper(params)
        predict_counts_with_params(
            full_params;
            cluster_model=cluster_model,
            emission_model=emission_model,
            redshift=redshift,
            predicted_bg_over_obs_time=predicted_obs_bg,
            shape=shape,
            pixel_edge_angle=pixel_edge_angle,
            observation_exposure_time=obs_exposure_time,
            response_function=response_function,
            centre_radius=centre_radius,
            mask=mask,
            integration_limit=integration_limit
        )
    end

    # a wrapper to handle running the gas model and likelihood calculation
    function likelihood_wrapper(params::AbstractVector{Float64})::Float64
        try
            return log_likelihood(
                observed,
                observed_background,
                predict_counts(params),
                predicted_bg_bg,
                log_obs_factorial
            )
        catch e
            if e isa PriorError || e isa ObservationError
                @mpidebug "Prior or observation error" e params
                return e.likelihood
            end
            rethrow()
        end
    end

    if plateau_tests > 0
        plateau_test(transform, prior_names, likelihood_wrapper; n_tests=plateau_tests)
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
            nsteps=length(prior_names),
            generate_direction=stepsampler.generate_mixture_random_direction,
            adaptive_nsteps="move-distance",
            max_nsteps=3 * length(prior_names),
            region_filter=true,
        )
    end

    # run Ultranest
    @mpiinfo "Launching sampler"
    ultranest_default_run_args = (
        max_num_improvement_loops=10,
        region_class=mlfriends.MLFriends,
    )
    @mpidebug "Ultranest arguments" ultranest_default_run_args ultranest_run_args
    merged = merge(ultranest_default_run_args, ultranest_run_args)
    results = sampler.run(;
        merged...
    )

    # output data
    @mpidebug "Sampler done"
    # print("result has these keys:", keys(results), "\n")
    sampler.print_results()
    sampler.plot_corner()
    sampler.plot_trace()
    sampler.plot_run()

    best_fit = results["maximum_likelihood"]["point"]
    best_fit_observation = predict_counts(best_fit)

    if MPI.Comm_rank(comm) == 0 && sampler.log == true
        try
            output_dir = sampler.logs["run_dir"]
            if output_dir isa AbstractString
                @mpiinfo "Running blob finder on best fit likelihood"
                p = run_blob_analysis(
                    observed,
                    log_likelihood_array(
                        observed,
                        observed_background,
                        best_fit_observation,
                        predicted_bg_bg,
                        log_obs_factorial
                    ),
                    get_centre_indices(
                        best_fit[1] * 1u"arcsecondᵃ",
                        best_fit[2] * 1u"arcsecondᵃ",
                        pixel_edge_angle,
                        redshift,
                        tuple(shape...)
                    ),
                )
                save("$output_dir/plots/blobs.svg", p)
            end
        catch e
            if !(e isa KeyError)
                rethrow(e)
            end
            @mpidebug "Skipping blob plot because unable to find log information"
        end
    end

    return (sampler, results, best_fit_observation)
end

"""
    sample(
        data::Dataset,
        energy_limits::NTuple{2, Unitful.Energy},
        cluster_model::Function,
        priors::AbstractVector{Prior},
        nhCol::SurfaceDensity,
        redshift::Real,
        x::NTuple{2,<:Real},
        y::NTuple{2,<:Real};
        bin_size::Real=10,
        use_interpolation::Bool=false,
        centre_radius=0,
        mask=nothing,
        cache_size::Int64=1000000000
        )

Run Bayesian inference on a given set of `data` considering only the selected
energy range.

* The cluster model can be any function that takes a set of parameters matching the priors and a `z` keyword argument, 
and returns two functions for the gas temperature and gas mass density as a function of radius.
* The first two priors should always be `x0` and `y0`, giving cluster centre position.
* `x` and `y` are tuples of `(min, max)` in pixels. These crop the observation.
* `mask` is optional. If included it should be a string pointing to a mask file using CIAO syntax. Only ellipses are supported.
* `centre_radius` excludes some radius, in pixels, around the centre from analysis
* `cache_size` controls the MEKAL cache size, in bytes (default is 1GB) [disabled]
* Additional kwargs will be passed through to the next `sample` function.
"""
function sample(
    data::Dataset,
    energy_limits::NTuple{2,<:Unitful.Energy},
    cluster_model::Function,
    priors::AbstractVector{<:Prior},
    nHcol::SurfaceDensity,
    redshift::Real,
    x::NTuple{2,<:Real},
    y::NTuple{2,<:Real};
    bin_size::Real=10,
    use_interpolation::Bool=false,
    centre_radius=0,
    mask=nothing,
    abundances=ones(10),
    kwargs...
)
    @argcheck [p.name for p in priors[1:2]] == ["x0", "y0"] || [p.name for p in priors[1:2]] == ["x", "y"]

    @mpiinfo "Loading response function"
    response_function, energy_range, channel_range = load_response(data, energy_limits...)

    # Check the number of energy bins in the response function matches the number of bins in the energy range
    @assert size(response_function, 2) == (length(energy_range) - 1) # subtract 1 because the range is bin edges

    @mpiinfo "Loading data"
    observation, observed_background = load_data(data)

    x_edges = x[1]:bin_size:x[2]
    y_edges = y[1]:bin_size:y[2]
    @mpidebug "Spatial bounds" x_edges y_edges

    @mpiinfo "Binning observation data"
    obs = bin_events(data, observation.first, channel_range, x_edges, y_edges)

    @mpiinfo "Binning background data"
    bg = bin_events(data, observed_background.first, channel_range, x_edges, y_edges)

    pixel_edge_angle = bin_size * data.pixel_edge_angle

    @mpidebug "Making transform"
    prior_names = [p.name for p in priors if !isa(p, DeltaPrior)]
    transform, param_wrapper = make_cube_transform(priors...)

    @mpiinfo "Generating emissions model"
    if use_interpolation
        emission_model = prepare_model_mekal_interpolation(nHcol, energy_range, redshift, abundances)
    else
        emission_model = prepare_model_mekal(nHcol, energy_range, redshift, abundances)
    end

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