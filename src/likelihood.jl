export DeltaPrior, LogUniformPrior, UniformPrior

"""
    log_likelihood(observed, observed_background, predicted, predicted_background, observed_log_factorial)

Calculate the log-likelihood of the prediction given an observation.

The observed and predicted arrays include background events.
Log factorial is calculated as `ln(observed) + ln(observed_background)`.
We require it to be supplied to improve performance - no need to calculate it every time.

We assume the predicted_background is scaled to the same exposure time as the observed background.
"""
function log_likelihood(
    observed,
    observed_background,
    predicted,
    predicted_background,
    observed_log_factorial
)
    @mpirankeddebug "Calculating log likelihood"

    @assert size(observed) == size(predicted) "Observations have size $(size(observed)) whereas predictions have size $(size(predicted))"
    @assert size(observed) == size(observed_background)
    @assert size(predicted) == size(predicted_background) || size(predicted_background) == ()


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

"""A delta prior that always returns a constant value."""
struct DeltaPrior{T<:Number} <: Prior
    value::T
end
function transform(prior::DeltaPrior, x::Real)
    @argcheck 0 <= x <= 1
    return prior.value
end

"""A uniform prior that draws from a uniform distribution between `min` and `max`."""
struct UniformPrior{T<:Number} <: Prior
    min::T
    max::T
    UniformPrior(min::T, max::T) where {T<:Number} = max > min ? new{T}(min, max) : error("Maximum is not greater than min")
end
function transform(prior::UniformPrior, x::Real)
    return x * (prior.max - prior.min) + prior.min
end

"""A log uniform prior that draws from a distribution between `min`
and `max` whose base 10 logarithm is uniformly distributed."""
struct LogUniformPrior{T<:Number} <: Prior
    min::T
    max::T
    LogUniformPrior(min::T, max::T) where {T<:Number} = max > min ? new{T}(min, max) : error("Maximum is not greater than min")
end
function transform(prior::LogUniformPrior, x::Real)
    lmin = log10(prior.min)
    lmax = log10(prior.max)
    return 10^(x * (lmax - lmin) + lmin)
end

"""
    make_cube_transform(priors::Prior...)

Turn a sequence of prior transform functions intoNorma a transform function that operates
on the hypercube generated by multinest.
"""
function make_cube_transform(priors::Prior...)

    """
        transform_cube(cube)

    Transforms the hypercube used by ultranest into physical prior values.

    Ultranest models priors as a unit hypercube where each dimesion is a unit uniform
    distribution. The transform function converts values on these uniform distributions
    to values on the physical prior distribution. Each column is a specific prior, so each
    row is a complete sample of the set of priors.
    """
    function transform_cube(cube::AbstractVector)
        # MT_200::Unitful.Mass,
        # fg_200,
        # a_GNFW,
        # b_GNFW,
        # c_GNFW,
        # c_500_GNFW,

        @mpirankeddebug "Transform started"

        for (c, p) in zip(axes(cube, 1), priors)
            cube[c] = transform.(Ref(p), cube[c])
        end

        @mpirankeddebug "Transform done"

        return cube
    end

    return transform_cube
end


