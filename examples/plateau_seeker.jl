using Random, PoissonRandom
using Unitful, UnitfulAstro, DimensionfulAngles
using BayesJ
using CairoMakie
using LinearAlgebra: I
using ProgressMeter
using DataFrames, PairPlots

z = 0.1

energy_bins = range(0.7u"keV", 6.0u"keV", step=0.01u"keV")
response_function = 250u"cm^2" * Matrix(I, length(energy_bins) - 1, length(energy_bins) - 1)

emission_model = BayesJ.prepare_model_mekal(
    0.0e20u"cm^-2",
    energy_bins,
    z,
)

function observe(temperature, density)
    let z = z, emission_model = emission_model, response_function = response_function
        return BayesJ.make_observation(
            temperature,
            density,
            z,
            (32, 32),
            0.492u"arcsecondᵃ",
            emission_model,
            3.0e6u"s",
            response_function,
            (0u"arcsecondᵃ", 0u"arcsecondᵃ"),
            0,
        )

    end
end

priors = [
    DeltaPrior("x0", 0.0),
    DeltaPrior("y0", 0.0),
    UniformPrior("n0", 0.1e-3, 40.0e-3),
    UniformPrior("n02", 0.001e-1, 6.0e-1),
    UniformPrior("rc", 100.0, 300.0),
    UniformPrior("rc2", 50.0, 200.0),
    UniformPrior("α", 0.1, 3.0),
    UniformPrior("β", 0.1, 2.0),
    UniformPrior("β2", 0.1, 5.0),
    UniformPrior("ϵ", 1.0, 5.0), # constrain ϵ<5
    UniformPrior("rs", 500.0, 1000.0),
    UniformPrior("T0", 1.0, 10.0),
    UniformPrior("Tmin/T0", 0.1, 1.0),
    UniformPrior("rcool", 100.0, 250.0),
    UniformPrior("acool", 1.0, 12.0),
    UniformPrior("rt", 1.0, 5.0),
    UniformPrior("a", 0.0, 1.0),
    UniformPrior("b", 1.0, 6.0),
    UniformPrior("c", 1.0, 4.0),
    DeltaPrior("d", 0.0),
]

prior_transform, param_wrapper = BayesJ.make_cube_transform(priors...)
prior_names = [p.name for p in priors if !isa(p, DeltaPrior)]
n_priors = length(prior_names)

n_tests = 100
max_tests = 1500
n = 0
t = 0

values = Array{Float64,2}(undef, n_priors, n_tests)

p = Progress(n_tests; dt=1.0)
while n < n_tests && t < max_tests
    global n
    global t

    t += 1

    rand_point = prior_transform(rand(n_priors))

    try
        model = Model_Vikhlinin2006(param_wrapper(rand_point)[3:end]..., z=z)
        observation = observe(model[1], model[2])
    catch e
        if e isa BayesJ.PriorError
            continue
        elseif e isa BayesJ.ObservationError
            n += 1
            values[:, n] = rand_point
            next!(p)
        else
            rethrow(e)
        end
    end
end

df = DataFrame(values[:, 1:n]', prior_names)

@info "Evaluations" n t

save("plateaus.svg", pairplot(df))