using Unitful, DimensionfulAngles
using Plots
using Profile

using PyCall
ultranest = pyimport_conda("ultranest", "ultranest", "conda-forge")
stepsampler = pyimport_conda("ultranest.stepsampler", "ultranest")

export sample

include("mpi.jl")
include("likelihood.jl")
include("gas_models.jl")
include("io.jl")

"""
    sample(observed::Array, observed_background::Array, response_function::Matrix, priors::PriorSet, obs_exposure_time, bg_exposure_time, redshift; emission_model, pixel_edge_angle, background_rate, average_effective_area, centre_radius)

Configure some necessary variables and launch ultranest.

* The observed array includes the background.
* The response function includes both the RMF and ARF, as described in `apply_response_function`.
* The emission model should be a function compatible with the requirements of the `surface_brightness` function, which it will be passed to.
* Prior names should match arguments for the gas model.
* The pixel edge angle describes the angular size observed by a single pixel in units such as arcseconds.
 This area is assumed to be square with the edge angle giving the side length.
* The average effective area is the effective area of the telescope averaged across energies,
 used with the total background rate across all channels (counts per unit telescope area per sky angle per second)
 to calculate the background counts per second per channel per pixel.
 * `centre_radius` controls how many pixels are excluded from the core.
"""
function sample(
    observed::T,
    observed_background::T,
    response_function::Matrix,
    priors::PriorSet,
    obs_exposure_time::Unitful.Time,
    bg_exposure_time::Unitful.Time,
    redshift::Real;
    emission_model,
    pixel_edge_angle=0.492u"arcsecondᵃ",
    background_rate=8.4e-6u"cm^-2/arcminuteᵃ^2/s",
    average_effective_area=250u"cm^2",
    centre_radius
) where {T<:AbstractArray}
    @mpidebug "Preparing for ultranest"

    predicted_bg_rate = background_rate / size(observed)[1] * average_effective_area * pixel_edge_angle^2
    predicted_obs_bg = predicted_bg_rate * obs_exposure_time # Used for adding background to observations
    predicted_bg_bg = predicted_bg_rate * bg_exposure_time # Used for log likelihood

    log_obs_factorial = log_factorial.(observed) + log_factorial.(observed_background)

    @assert all(isfinite, observed)
    @assert all(isfinite, observed_background)
    @assert all(isfinite, log_obs_factorial)

    @assert size(observed) == size(observed_background)

    shape = [i for i in size(observed)][2:3]

    @mpiinfo "Observation has shape $(size(observed))"
    @mpiinfo "Background has shape $(size(observed_background))"
    @mpiinfo "Response matrix has shape $(size(response_function))"

    # a wrapper to handle running the gas model and likelihood calculation
    @mpidebug "Generating likelihood wrapper"
    function likelihood_wrapper(params)
        @mpirankeddebug "Likelihood wrapper called" params

        prior_pairs = priors.cube_to_name(params)::Dict{<:AbstractString,<:Prior}

        predicted = Model_NFW_GNFW(;
            prior_pairs...,
            redshift=redshift,
            shape=shape,
            pixel_edge_angle=pixel_edge_angle,
            emission_model=emission_model,
            obs_exposure_time=obs_exposure_time,
            response_function=response_function,
            centre_radius=centre_radius
        )

        predicted .+= predicted_obs_bg

        @mpirankeddebug "Predicted results generated"

        return log_likelihood(
            observed,
            observed_background,
            predicted,
            predicted_bg_bg,
            log_obs_factorial
        )
    end

    # ultranest setup
    @mpidebug "Creating sampler"
    # paramnames = ["MT_200", "fg_200", "x_0", "y_0"] # move to pairs with prior objects?
    paramnames = ["MT_200", "fg_200"]
    sampler = ultranest.ReactiveNestedSampler(
        paramnames,
        likelihood_wrapper,
        transform=priors.transform,
        vectorized=false,
        log_dir="logs"
    )

    # @mpidebug "Creating stepsampler"
    # sampler.stepsampler = stepsampler.SliceSampler(
    #     nsteps=25,
    #     generate_direction=stepsampler.generate_mixture_random_direction,
    # )

    # run Ultranest
    @mpiinfo "Launching sampler"
    results = sampler.run()

    # output data
    @mpidebug "Sampler done"
    # print("result has these keys:", keys(results), "\n")
    sampler.print_results()
    sampler.plot_corner()
    sampler.plot_trace()
    sampler.plot_run()

    return (sampler, results)
end

"""
    sample(data::Dataset, energy_range, priors::Dict, nhCol, redshift; bin_size, use_interpolation, centre_radius)

Run Bayesian inference on a given set of `data`, considering only the selected energy range.

* `bin_size` controls spatial size of bins, in pixels
* `use_interpolation` controls whether the gas emission model uses interpolation
* `centre_radius` controls the excluded area in the core, in pixels
"""
function sample(
    data::Dataset,
    energy_range::AbstractRange{<:Unitful.Energy},
    priors,
    nHcol::SurfaceDensity,
    redshift::Real;
    bin_size::Real=10,
    use_interpolation::Bool=true,
    centre_radius=4
)
    @mpiinfo "Loading data"

    observation, observed_background = load_data(data)

    obs = bin_events(data, observation.first, energy_range, 3700:bin_size:4200, 4100:bin_size:4550)
    bg = bin_events(data, observed_background.first, energy_range, 3700:bin_size:4200, 4100:bin_size:4550)
    pixel_edge_angle = bin_size * data.pixel_edge_angle
    @mpidebug "Done binning events"

    @assert size(obs) == size(bg)

    @mpidebug "Making transform"
    prior_set = generate_transform(priors)

    @mpidebug "Calling load_response"
    response_function = load_response(data, energy_range)

    @mpiinfo "Generating emissions model"
    emission_model = prepare_model_mekal(nHcol, energy_range, redshift, use_interpolation=use_interpolation)

    @mpiinfo "Testing emissions model"
    em_direct = prepare_model_mekal(nHcol, energy_range, redshift, use_interpolation=false)

    model = Model_NFW_GNFW(
        MT_200=5e14u"Msun",
        fg_200=0.13,
        α=1.0510, # Using universal values from Arnaud 2010
        β=5.4905,
        γ=0.3081,
        c_500_GNFW=1.177,
        z=redshift,
        shape=[64, 64],
        pixel_edge_angle=0.492u"arcsecondᵃ",
        emission_model=emission_model,
        exposure_time=100e3u"s",
        response_function=response_function,
        centre_coordinates=(0u"arcsecondᵃ", 0u"arcsecondᵃ"),
        centre_radius=centre_radius
    )
    model_direct = model = Model_NFW_GNFW(
        MT_200=5e14u"Msun",
        fg_200=0.13,
        α=1.0510, # Using universal values from Arnaud 2010
        β=5.4905,
        γ=0.3081,
        c_500_GNFW=1.177,
        z=redshift,
        shape=[64, 64],
        pixel_edge_angle=0.492u"arcsecondᵃ",
        emission_model=em_direct,
        exposure_time=100e3u"s",
        response_function=response_function,
        centre_coordinates=(0u"arcsecondᵃ", 0u"arcsecondᵃ"),
        centre_radius=centre_radius
    )
    replace!(model, NaN => 0)
    replace!(model_direct, NaN => 0)

    err = sum(abs2, model - model_direct)
    @mpiinfo "Error in emission model is" err

    sample(
        obs,
        bg,
        response_function,
        prior_set,
        observation.second,
        observed_background.second,
        redshift;
        emission_model=emission_model,
        pixel_edge_angle=pixel_edge_angle,
        centre_radius=centre_radius
    )
end