using PyCall

ultranest = pyimport_conda("ultranest", "ultranest")

include("gas_models.jl")

function log_likelihood(observed::Array{T}, predicted::Array{N}) where {T<:Real,N<:Real}
    t1 = observed .* log.(predicted) - predicted - logCobs_factorial
    # TODO: Background

    replace!(t1, -Inf => 0)

    return sum(t1)
end

function transform(cube)
    # MT_200::Unitful.Mass,
    # fg_200,
    # a_GNFW,
    # b_GNFW,
    # c_GNFW,
    # c_500_GNFW,

    @debug "Transform started"

    # MT_200
    cube[1] = cube[1] * (1e15 - 1e14) + 1e14
    # fg_200
    cube[2] = cube[2] * (0.2 - 0.05) + 0.05
    # a_GNFW
    cube[3] = 1.062
    # b_GNFW
    cube[4] = 5.4807
    # c_GNFW
    cube[5] = 0.3292
    # c_500_GNFW
    cube[6] = 1.156

    @debug "Transform done"

    return cube
end

const dshape = [12, 12]
const observed = ceil.(Int64, complete_matrix(Model_NFW_GNFW(
        5e14u"Msun",
        0.13,
        1.0620,
        5.4807,
        0.3292,
        1.156,
        0.1,
        dshape,
        0.492u"arcsecond"
    ), dshape
))
log_fact(n::Integer) = sum(log.(1:n))
logCobs_factorial = log_fact.(observed)

function likelihood_wrapper(params)
    @debug "Likelihood started"
    n, d = size(params)
    params_rows = Vector{Vector{eltype(params)}}(undef, n)
    for i in 1:n
        params_rows[i] = params[i, :]
    end
    # display(params_rows)
    predicted = complete_matrix.([Model_NFW_GNFW(params_rows[i]...,
            0.1,
            dshape,
            0.492u"arcsecond"
        ) for i in 1:n], Ref(dshape))
    @debug "Predicted results generated"
    return log_likelihood.(Ref(observed), predicted)
end

function run_ultranest()
    paramnames = ["MT_200", "fg_200", "a", "b", "c", "c_500_GNFW"]
    sampler = ultranest.ReactiveNestedSampler(paramnames, likelihood_wrapper, transform=transform, vectorized=true)
    @debug "Sampler starting"
    results = sampler.run()
    @debug "Sampler done"
    print("result has these keys:", keys(results), "\n")
    sampler.print_results()
    sampler.plot()
end

run_ultranest()