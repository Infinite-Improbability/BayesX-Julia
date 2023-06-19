using PyCall
ultranest = pyimport_conda("ultranest", "ultranest", "conda-forge")

include("gas_models.jl")

"""
    log_likelihood(observed, observed_background, predicted, predicted_background, observed_log_factorial)

Calculate the log-likelihood of the prediction given an observation.

The observed and predicted arrays include background events.
Log factorial is calculated as `ln(observed) + ln(observed_background)`.
We require it to be supplied to improve performance - no need to calculate it every time.

"""
function log_likelihood(
    observed::A{N},
    observed_background::A{N},
    predicted::A{T},
    predicted_background::A{T},
    observed_log_factorial
) where {A<:AbstractArray,T<:Real,N<:Integer}
    @. t1 = observed * log(predicted) - predicted +
            observed_background .* log(predicted_background) - predicted_background -
            observed_log_factorial

    # If we have zero counts we can't take the log so
    # we'll just skip over it.
    replace!(t1, -Inf => 0)

    return sum(t1)
end

"""
    transform(cube)

Transforms the hypercube used by ultranest into physical prior values.

Ultranest models priors as a unit hypercube where each dimesion is a unit uniform
distribution. The transform function converts values on these uniform distributions
to values on the physical prior distribution. Each column is a specific prior, so each
row is a complete sample of the set of priors.
"""
function transform(cube::A) where {A<:AbstractArray}
    # MT_200::Unitful.Mass,
    # fg_200,
    # a_GNFW,
    # b_GNFW,
    # c_GNFW,
    # c_500_GNFW,

    @debug "Transform started"

    # MT_200
    @. cube[:, 1] = cube[:, 1] * (1e15 - 1e14) + 1e14
    # fg_200
    @. cube[:, 2] = cube[:, 2] * (0.2 - 0.05) + 0.05

    @debug "Transform done"

    return cube
end

const surrogate_model = prepare_model(2.2, 0.1, 0.3:0.1:3.0)
const dshape = [24, 24]
const observed = round.(Int64, complete_matrix(Model_NFW_GNFW(
        5e14u"Msun",
        0.13,
        1.0620,
        5.4807,
        0.3292,
        1.156,
        0.1,
        dshape,
        0.492u"arcsecond",
        surrogate_model
    ), dshape
))

"""
    log_factorial(n)

    Finds the natural logarithm of the factorial of `n`.

    `n` rapidly gets to large to quickly and directly calculate the factorial
    so we exploit logarithm rules to expand it out to a series of sums.
"""
log_factorial(n::N) where {N<:Integer} = sum(log.(1:n))

const logCobs_factorial = log_factorial.(observed)

function likelihood_wrapper(params)
    @debug "Likelihood started"
    n, _ = size(params)
    params_rows = Vector{Vector{eltype(params)}}(undef, n)
    for i in 1:n
        params_rows[i] = params[i, :]
    end
    # display(params_rows)
    predicted = complete_matrix.(
        [
            Model_NFW_GNFW(params_rows[i]...,
                1.062,
                5.4807,
                0.3292,
                1.156,
                0.1,
                dshape,
                0.492u"arcsecond",
                surrogate_model,
            ) for i in 1:n
        ],
        Ref(dshape)
    )
    @debug "Predicted results generated"
    return log_likelihood.(Ref(observed), predicted)
end

function run_ultranest()
    paramnames = ["MT_200", "fg_200"]
    sampler = ultranest.ReactiveNestedSampler(paramnames, likelihood_wrapper, transform=transform, vectorized=true)
    @debug "Sampler starting"
    results = sampler.run()
    @debug "Sampler done"
    print("result has these keys:", keys(results), "\n")
    sampler.print_results()
    sampler.plot()
    return (sampler, results)
end

s, r = run_ultranest()