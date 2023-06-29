include("bayes.jl")

abstract type BayesXDataset end
# Include exposure times in dataset
# And NHcol, bg count rate?

"""
    run(data, energy_range)

Run Bayesian inference on a given set of `data`, considering only the selected
energy range. An gas emission model `(density, temperature) → emissivity` can be provided.
"""
function run(
    data::BayesXDataset,
    energy_range,
    priors,
    ;
    emission_model = prepare_model_mekal(2.2, 0.1, 0.3:0.1:3.0),
    pixel_size_length = 0.492u"arcsecond"
    )