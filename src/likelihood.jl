using Distributions
using SpecialFunctions: loggamma

export DeltaPrior, LogUniformPrior, UniformPrior, NormalPrior, GenericPrior, DependentUniformPrior, DependentLogUniformPrior

"""
    log_likelihood(observed, observed_background, predicted, predicted_background, observed_log_factorial)

Calculate the log-likelihood of the prediction given an observation.

The observed and predicted arrays include background events.

The constant likelihood includes all parts of the likelihood that remain
constant across observations. We require it to be supplied to improve performance - no need to calculate it every time. 
It also includes the log factorial components calculated as `ln(observed) + ln(observed_background)`.
As we haven't implemented support for background fitting on the user input side of things this includes the background.
"""
function log_likelihood(
    observed,
    predicted,
    constant_likelihood
)
    # @mpirankeddebug "Calculating log likelihood"

    t = log_likelihood_array(
        observed,
        predicted,
        constant_likelihood
    )
    # @mpirankeddebug "likelihood is" sum(skipmissing(t))
    return sum(skipmissing(t))
end

function log_likelihood_array(
    observed,
    predicted,
    constant_likelihood
)
    t = @. observed * log(predicted) - predicted + constant_likelihood
    @assert all(isfinite, skipmissing(t))

    return t
end

"""
    log_factorial(n)

Finds ``\\ln{n!}`` the natural logarithm of the factorial of `n`.

As Γ(n) = (n-1)! we can use SpecialFunctions.jl loggamma function to quickly calculate the log factorial.

It is intended to be broadcast across all values of the data array.
"""
log_factorial(n::N) where {N<:Number} = loggamma(n + 1)

"""
Abstract supertype for priors. Should implement:
- A transform(prior, x) function that transforms a value x on the unit range to a value on the distribution represented by the prior.
- A field `name<:AbstractString` 
"""
abstract type Prior end

"""
    transform(prior, x)

Transforms a value x on the unit range to a value on the distribution represented by the prior.
"""
function transform(prior::Prior, x::Real) end

"""
    DeltaPrior(name::AbstractString, value::Number)

A delta prior that always returns a constant value.
"""
struct DeltaPrior{S<:AbstractString,T<:Number} <: Prior
    name::S
    value::T
end
function transform(prior::DeltaPrior, ::Real)
    return prior.value
end

"""
    UniformPrior(name::AbstractString, min::Number, max::Number)

A uniform prior that draws from a uniform distribution between `min` and `max`.
"""
struct UniformPrior{S<:AbstractString,T<:Number} <: Prior
    name::S
    min::T
    max::T
    UniformPrior(name::S, min::T, max::T) where {S<:AbstractString,T<:Number} = max > min ? new{S,T}(name, min, max) : error("Maximum is not greater than min for prior $name")
end
function transform(prior::UniformPrior, x::Real)
    return x * (prior.max - prior.min) + prior.min
end

"""
    LogUniformPrior(name::AbstractString, min::Number, max::Number)

A log uniform prior that draws from a distribution between `min` and `max` whose base 10 logarithm is uniformly distributed."""
struct LogUniformPrior{S<:AbstractString,T<:Number} <: Prior
    name::S
    min::T
    max::T
    LogUniformPrior(name::S, min::T, max::T) where {S<:AbstractString,T<:Number} = max > min ? new{S,T}(name, min, max) : error("Maximum is not greater than min for prior $name")
end
function transform(prior::LogUniformPrior, x::Real)
    lmin = log10(prior.min)
    lmax = log10(prior.max)
    return 10^(x * (lmax - lmin) + lmin)
end

"""
    NormalPrior(name::AbstractString, mean::Number, σ::Number)

A normal prior that draws from a Gaussian/normal distribution with `mean` and standard deviation `σ`.
"""
struct NormalPrior{S<:AbstractString,T<:Number} <: Prior
    name::S
    mean::T
    σ::T
    dist::Normal{T}
end
NormalPrior(name::S, mean::T, σ::T) where {S<:AbstractString,T<:Number} = NormalPrior(name, mean, σ, Normal(mean, σ))
function transform(prior::NormalPrior, x::Real)
    return quantile(prior.dist, x)
end

"""
    GenericPrior(name::AbstractString, dist::Distributions.UnivariateDistribution)

A generic prior that draws from the specified distribution.

See [Distributions.jl](https://juliastats.org/Distributions.jl/stable/univariate/) for a list of avaliable distributions.
"""
struct GenericPrior{S<:AbstractString,D<:UnivariateDistribution} <: Prior
    name::S
    dist::UnivariateDistribution
end
function transform(prior::GenericPrior, x::Real)
    return quantile(prior.dist, x)
end

abstract type DependentPrior <: Prior end

function transform(::DependentPrior, x::Real)
    return x
end

"""
    DependentUniformPrior(name::AbstractString, depends_on::AbstractString, range::Number)

Generates a dependent uniform prior.

The dependent prior has a minimum value equal to that of the parent prior and a maximum value equal to the parent prior plus the range.
By specifying a negative range, the dependent prior can be made to have a maximum value less than the parent prior.
"""
struct DependentUniformPrior{S<:AbstractString,T<:Number} <: DependentPrior
    name::S
    depends_on::S
    range::T
end
function transform(prior::DependentUniformPrior, x::Real, parent_value::Real)
    return x * prior.range + parent_value
end

"""
    DependentLogUniformPrior(name::AbstractString, depends_on::AbstractString, range::Number)

Generates a dependent log uniform prior.

The dependent prior has a minimum value equal to that of the parent prior and a maximum value equal to the parent prior plus the range.
By specifying a negative range, the dependent prior can be made to have a maximum value less than the parent prior.
"""
struct DependentLogUniformPrior{S<:AbstractString,T<:Number} <: DependentPrior
    name::S
    depends_on::S
    range::T
end
function transform(prior::DependentLogUniformPrior, x::Real, parent_value::Real)
    max = parent_value + prior.range
    lmin = log10(parent_value)
    lmax = log10(max)

    return 10^(x * (lmax - lmin) + lmin)
end


"""
    make_cube_transform(priors::Prior...)::NTuple{2, Function}

Turn a sequence of prior transform functions into a transform function `transform_cube` that operates
on the hypercube generated by multinest, stripping delta functions.

It also returns a function `reconstruct_args` take takes the output of the first function and returns it with
the values of the delta priors inserted in appropriate postions. 
This is because Ultranest can run with delta priors but outputs error warnings in the process.

Returns `transform_cube` and `reconstruct_args` as a tuple.
"""
function make_cube_transform(priors::Prior...)::NTuple{2,Function}

    delta_priors = DeltaPrior[]
    variable_priors = Prior[]
    is_delta = Bool[]
    dependencies = []

    for p in priors
        if isa(p, DeltaPrior)
            push!(delta_priors, p)
            push!(is_delta, true)
        else
            push!(variable_priors, p)
            push!(is_delta, false)
        end
    end

    for (i, p) in enumerate(variable_priors)
        if isa(p, DependentPrior)
            depends_on = findfirst(x -> x.name == p.depends_on, variable_priors)
            if isnothing(depends_on)
                @mpierror "Prior $(p.name) depends on $(p.depends_on) but no such non-delta prior exists"
                throw(ErrorException("Unable to generate prior transform"))
            end
            @mpidebug "Prior $(p.name) ($i) depends on $(p.depends_on) ($depends_on)"
            push!(dependencies, (i, depends_on))
        end
    end

    @assert length(is_delta) == length(variable_priors) + length(delta_priors)
    @assert count(is_delta) == length(delta_priors)
    @assert count(is_delta .== false) == length(variable_priors)

    """
        transform_cube(cube)

    Transforms the hypercube used by ultranest into physical prior values.

    Ultranest models priors as a unit hypercube where each dimesion is a unit uniform
    distribution. The transform function converts values on these uniform distributions
    to values on the physical prior distribution. Each column is a specific prior, so each
    row is a complete sample of the set of priors.
    """
    function transform_cube(cube::AbstractVector)

        @mpirankeddebug "Transforming prior hypercube"

        for (c, p) in zip(axes(cube, 1), variable_priors)
            cube[c] = transform.(Ref(p), cube[c])
        end

        return cube
    end

    transform_wrapper = transform_cube
    if length(dependencies) > 0
        function transform_dependents(cube::AbstractVector)
            cube = transform_cube(cube)
            for (prior, depends_on) in dependencies
                cube[prior] = transform(variable_priors[prior], cube[prior], cube[depends_on])
            end
            return cube
        end
        transform_wrapper = transform_dependents
    end

    delta_values = [d.value for d in delta_priors]

    """
        reconstruct_args(unfixed::AbstractVector)

    Takes a list of values generated by transforming the unit hypercube and inserts values associated with
     the delta priors in the appropriate places.

    This is required because Ultranest does not handle delta priors well.
    """
    function reconstruct_args(unfixed::AbstractVector)
        let dv = delta_values, id = is_delta

            dv_i = 0 # delta prior index
            un_i = 0 # unfixed prior index

            Tuple(
                i ? begin # if delta insert next value from delta values
                    dv_i += 1
                    dv[dv_i]
                end : begin
                    un_i += 1 # else insert next unfixed prior
                    unfixed[un_i]
                end
                for i in id
            )
        end
    end

    # Verify that our reconstructed results match the expected values
    @assert reconstruct_args([transform(p, 0.5) for p in variable_priors]) == Tuple(transform(p, 0.5) for p in priors)

    return transform_wrapper, reconstruct_args
end


