using PyCall
ultranest = pyimport_conda("ultranest", "ultranest", "conda-forge")

include("gas_models.jl")
include("io.jl")

"""
    log_likelihood(observed, observed_background, predicted, predicted_background, observed_log_factorial)

Calculate the log-likelihood of the prediction given an observation.

The observed and predicted arrays include background events.
Log factorial is calculated as `ln(observed) + ln(observed_background)`.
We require it to be supplied to improve performance - no need to calculate it every time.
"""
function log_likelihood(
    observed,
    observed_background,
    predicted,
    predicted_background,
    observed_log_factorial
)
    @assert size(observed) == size(predicted) "Observations have size $(size(observed)) whereas predictions have size $(size(predicted))"
    @assert size(observed) == size(observed_background)
    @assert size(predicted) == size(predicted_background)


    t1 = @. observed * log(predicted) - predicted
    t2 = @. observed_background * log(predicted_background) - predicted_background

    replace!(t2, NaN => 0) # for some reason log(0) sometimes produces NaN not -Inf

    t = @. t1 + t2 - observed_log_factorial

    # If we have zero counts we can't take the log so
    # we'll just skip over it.
    replace!(t, -Inf => 0, NaN => 0)

    # For some reason 

    @assert all(isfinite, t) display(t)

    return sum(t)
end

"""
    log_factorial(n)

Finds the natural logarithm of the factorial of `n`.

`n` rapidly gets to large to quickly and directly calculate the factorial
so we exploit logarithm rules to expand it out to a series of sums.

It is intended to be broadcast across all values of the data array.
"""
log_factorial(n::N) where {N<:Integer} = sum(log.(1:n))

"""
Abstract supertype for priors. Should implement a transform(prior, x) function that Transforms
a value x on the unit range to a value on the distribution represented by the prior.
"""
abstract type Prior end

"""
    transform(prior, x)
Transforms a value x on the unit range to a value on the distribution represented by the prior.
"""
function transform(prior::Prior, x::Real) end

"""A delta prior, that always returns a constant value."""
struct DeltaPrior{T<:Number} <: Prior
    value::T
end
function transform(prior::DeltaPrior, x::Real)
    @argcheck 0 <= x <= 1
    return prior.value
end

"""A uniform prior, that draws from a uniform distribton between `min` and `max`."""
struct UniformPrior{T<:Number} <: Prior
    min::T
    max::T
    UniformPrior(min::T, max::T) where {T<:Number} = max > min ? new{T}(min, max) : error("Maximum is not greater than min")
end
function transform(prior::UniformPrior, x::Real)
    return x * (prior.max - prior.min) + prior.min
end

function make_cube_transform(priors::Prior...)

    """
        transform_cube(cube)

    Transforms the hypercube used by ultranest into physical prior values.

    Ultranest models priors as a unit hypercube where each dimesion is a unit uniform
    distribution. The transform function converts values on these uniform distributions
    to values on the physical prior distribution. Each column is a specific prior, so each
    row is a complete sample of the set of priors.
    """
    function transform_cube(cube::AbstractArray)
        # MT_200::Unitful.Mass,
        # fg_200,
        # a_GNFW,
        # b_GNFW,
        # c_GNFW,
        # c_500_GNFW,

        @debug "Transform started"

        for (c, p) in zip(axes(cube, 2), priors)
            cube[:, c] = transform.(Ref(p), cube[:, c])
        end

        @debug "Transform done"

        return cube
    end

    return transform_cube
end

"""
    run(observed, observed_background, response_function, transform, exposure_time; emission_model, pixel_edge_angle)

Configure some necessary variables and launch ultranest. The observed array includes
the background. The model is an interpolation over a true emission model.
"""
function run(
    observed::T,
    observed_background::T,
    response_function::Matrix,
    transform::Function,
    exposure_time::Unitful.Time;
    emission_model,
    pixel_edge_angle=0.492u"arcsecond"
) where {T<:AbstractArray}
    # TODO: actual background predictions

    predicted_bg = observed_background * 0 .+ 1

    log_obs_factorial = log_factorial.(observed) + log_factorial.(observed_background)

    @assert all(isfinite, observed)
    @assert all(isfinite, observed_background)
    @assert all(isfinite, log_obs_factorial)

    @assert size(observed) == size(observed_background)

    dshape = [i for i in size(observed)][2:3]

    # a wrapper to handle running the gas model and likelihood calculation
    function likelihood_wrapper(params)
        @debug "Likelihood started"

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
            0.1,
            dshape,
            pixel_edge_angle,
            emission_model,
            exposure_time,
            response_function
        ) for i in 1:n]


        @debug "Predicted results generated"

        return log_likelihood.(
            Ref(observed),
            Ref(observed_background),
            predicted,
            Ref(predicted_bg),
            Ref(log_obs_factorial)
        )
    end

    # ultranest setup
    paramnames = ["MT_200", "fg_200"]
    sampler = ultranest.ReactiveNestedSampler(
        paramnames,
        likelihood_wrapper,
        transform=transform,
        vectorized=true
    )

    # run Ultranest
    @debug "Sampler starting"
    results = sampler.run()
    @debug "Sampler done"

    # print("result has these keys:", keys(results), "\n")

    # output data
    sampler.print_results()
    # sampler.plot()

    return (sampler, results)
end

"""
    run(data::BayesXDataset, energy_range, priors)

Run Bayesian inference on a given set of `data`, considering only the selected
energy range. An gas emission model `(density, temperature) → emissivity` can be provided.
"""
function run(
    data::BayesXDataset,
    energy_range,
    priors;
    nHcol=2.2 # units of 10²² atoms per cm⁻²
)

    observation, observed_background = load_data(data)

    obs = bin_events(observation, energy_range, 2000:100:4000, 2000:100:4000)
    bg = bin_events(observed_background, energy_range, 2000:100:4000, 2000:100:4000)

    @assert size(obs) == size(bg)

    transform = make_cube_transform(priors...)

    response_function = load_response(data, energy_range)

    emission_model = prepare_model_mekal(nHcol, 0.1, LinRange(energy_range[1], energy_range[2], size(response_function)[2] + 1)) # we need this +1 but it seems to be one element too short

    run(obs, bg, response_function, transform, data.exposure_time; emission_model=emission_model, pixel_edge_angle=data.pixel_edge_angle)
end

data = FITSData("/home/ryan/data/chandra/4361/manual3/repro/acisf04361_repro_evt2.fits", "/home/ryan/data/chandra/4361/manual3/repro/bg_trimmed_300-7000.fits", "/home/ryan/data/chandra/4361/manual3/repro/specx/specx.arf", "/home/ryan/data/chandra/4361/manual3/repro/specx/specx.rmf", 300000u"s", 0.492u"arcsecond")
priors = [UniformPrior(1.0e14, 1.0e15), UniformPrior(0.08, 0.2)]
#nHcol = 3.89