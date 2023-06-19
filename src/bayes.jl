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
    observed,
    observed_background,
    predicted,
    predicted_background,
    observed_log_factorial
)
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

"""
    log_factorial(n)

Finds the natural logarithm of the factorial of `n`.

`n` rapidly gets to large to quickly and directly calculate the factorial
so we exploit logarithm rules to expand it out to a series of sums.

It is intended to be broadcast across all values of the data array.
"""
log_factorial(n::N) where {N<:Integer} = sum(log.(1:n))


"""
    run_ultranest(observed, observed_background, model)

Configure some necessary variables and launch ultranest. The observed array includes
the background. The model is an interpolation over a true emission model.
"""
function run_ultranest(
    observed::T,
    observed_background::T,
    model=prepare_model_mekal(2.2, 0.1, 0.3:0.1:3.0),
) where {T<:AbstractArray}

    # shape of the data
    # dshape = [i for i in size(observed)][1:2]

    # fake a background
    # TODO: actual background predictions
    predicted_bg = observed_background * 1.0

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
                    model,
                ) for i in 1:n
            ],
            Ref(dshape)
        )

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
    @info "Sampler starting"
    results = sampler.run()
    @info "Sampler done"

    # print("result has these keys:", keys(results), "\n")

    # output data
    sampler.print_results()
    sampler.plot()

    return (sampler, results)
end

const sh = [12, 12]
const m = prepare_model_mekal(2.2, 0.1, 0.3:0.1:3.0)
const obs = round.(Int64, complete_matrix(Model_NFW_GNFW(
        5e14u"Msun",
        0.13,
        1.0620,
        5.4807,
        0.3292,
        1.156,
        0.1,
        sh,
        0.492u"arcsecond",
        m
    ), sh
))

s, r = run_ultranest(
    obs,
    obs * 0,
);