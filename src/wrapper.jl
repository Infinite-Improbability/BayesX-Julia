include("bayes.jl")
include("io.jl")


"""
    run(data, energy_range)

Run Bayesian inference on a given set of `data`, considering only the selected
energy range. An gas emission model `(density, temperature) → emissivity` can be provided.
"""
function run(
    data::BayesXDataset,
    energy_range,
    priors;
    emission_model=prepare_model_mekal(2.2, 0.1, 0.3:0.1:3.0),
    pixel_size_length=0.492u"arcsecond"
)
    observation, observed_background = load_data(data) # TODO: binning

    bin_events(obs, energy_range, 2000:4000, 2000:4000)

    transform = make_cube_transform(priors)

    _run_ultranest(observation, observed_background, emission_model)


    # Apply energy range restrictions and bin data
    # Turn priors into transform function
    # Pass through to sampler

end